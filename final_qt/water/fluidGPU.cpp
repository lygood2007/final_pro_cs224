/**
 *@NOTE the following assumptions:
 *  all distances in meters
 *  all time steps in seconds
 *  gravity is in the -y plane
 *  2D height field plane is the xz plane
 *
 *Comment history:
 *  04/14 - added comments for the algo structure from the paper - SH
 *
 *
**/

#include "fluidGPU.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <QImage>
#include <QColor>
#include <QString>
#include <string>
#include <string.h>
#include <iostream>
#include <sstream>
#include "utils.h"
#include "CS123Common.h"
#include "debug_marco.h"

using std::cout;
using std::endl;

extern "C"
{
void initGridGPU( const int hostGridSize, const int hostGridPaintSize, const float hostdx, const float halfdm, const float* hostTerrainMap );
void copybackGPU(FieldType type, float* hostMap  );
void destroyGPUmem();
void addDropGPU(const int posX, const int posZ, const int radius, const float h );
void advectGPU(const float dt);
void updateFluidGPU( const float dt );
bool findSupportDevice();

void initParticlesGPU(const float minHeight, const int numSplashParticles, const int numSprayParticles, const int numFoamParticles);
void updateParticlesGPU( const float minHeight, const float dt, const float halfDomain, const float mdxInv, const float accX, const float accY, const float accZ );
void intersectParticlesGPU( const float minHeight, const float halfDomain, const float mdx, const float mdxInv,
                            const float splashVeff, const float splashHeightChange, const float sprayVeff, const float sprayHeightChange );
void inputParticlesGPU( const float *particlePositions, const float *particleVelocities );
void inputSprayParticlesGPU( const float *particlePositions, const float *particleVelocities );
void inputFoamParticlesGPU( const float *particlePositions, const float *ttlArray );

void checkBreakingWavesGPU( const float condition1, const float condition2, const float condition3,
                            const float mdx, const float mdxInv, const float dt, const float halfDomain,
                            const float heightChange, const float lambdaY, const int numParticlesToAdd );
void inputDepthGPU( const float* newDepthField );

void clampFieldsGPU( const float velocityClamp );

void initDampeningFieldsGPU( const int sizeDampeningRegion, const float quadraticA, const float quadraticB, const float quadraticC );
void dampenWavesGPU( const float hRest, const float dt, const float dxInv, const float lambdaUpdate, const float lambdaDecay );

void glBindBuffer (GLenum target, GLuint buffer);
void  glGenBuffers (GLsizei n, GLuint *buffers);
void *glMapBuffer(	GLenum target,GLenum access);
void  glBufferData (GLenum target, GLsizeiptr size, const GLvoid *data, GLenum usage);
void *glUnmapBuffer(GLenum target);
}

FluidGPU::FluidGPU()
{
    // Default initialization
    init();
}

FluidGPU::FluidGPU(const int gridSize, const float domainSize)
{
    init( gridSize, domainSize );
}

FluidGPU::FluidGPU(Terrain *t , GLWidget *glw)
{
    m_glw = glw;
    init( t->getGridLength(), 2*t->getBound() );

}

FluidGPU::~FluidGPU()
{
    // Release the heap
    safeFreeArray1D( m_paintField );
    safeFreeArray1D( m_indices );
    safeFreeArray1D( m_velocityU );
    safeFreeArray1D( m_velocityW );
    safeFreeArray1D( m_terrainHeightField );
    safeFreeArray1D( m_sigmaField );
    safeFreeArray1D( m_gammaField );
    safeFreeArray1D( m_phiField );
    safeFreeArray1D( m_psiField );
    safeFreeArray1D( m_paintNormalField );
    safeFreeArray1D( m_heightField );
    safeFreeArray1D( m_depthField );

    if(m_glw->m_useParticles)
    {
        safeFreeArray1D( m_spray_positions );
        safeFreeArray1D( m_spray_velocities );
        safeFreeArray1D( m_splash_positions );
        safeFreeArray1D( m_splash_velocities );
        safeFreeArray1D( m_foam_positions );
        safeFreeArray1D( m_foam_ttls );
        safeFreeArray1D( m_splash_to_foam_array );
    }

    destroyGPUmem();

    // release the particles
    for(int i = 0; i < m_particles.size(); i++){
        Particle *currParticle = m_particles[i];
        m_particles[i] = NULL;
        delete currParticle;
    }
    // release the particle sources
    for(int i = 0; i < m_particle_sources.size(); i++){
        ParticleSource *currParticleSource = m_particle_sources[i];
        m_particle_sources[i] = NULL;
        delete currParticleSource;
    }
}

//where the magic happens
void FluidGPU::draw() const
{
    drawFluid( DRAW_MESH_VBO );
#ifdef USE_PARTICLES
    //  drawParticles();
#endif

//    if(m_glw->m_useParticles) drawParticles2(); //moved this call to glwidget for rendering separatly

    if( m_renderNormals )
        drawNormal();
}

/**
 * Update the simulation at each time step
**/
void FluidGPU::update(const float dt)
{
    if( dt > 0.05 )
    {
        // We assume this is not stable
        // In fact, delta t should be smaller than delta_x/(g*D), where D is the maximum depth of the fluid
        return;
    }

    m_dt = dt;

    updateFluidGPU( dt );
   // clampFields();
    copybackGPU(HEIGHT,m_heightField);

    copybackGPU(VEL_U,m_velocityU);
    copybackGPU(VEL_W,m_velocityW);
    copybackGPU( PAINT, (float*) m_paintField );
    copybackGPU(NORMAL,(float*)m_paintNormalField);

#ifdef USE_PARTICLES
    //updateParticles();
    //  updateParticleSources();
#endif

    if(m_glw->m_useParticles)
    {
        float halfDomain = m_domainSize / 2.0;

        updateParticlesGPU( TERRAIN_MIN_HEIGHT, dt, halfDomain, m_dxInv, m_particle_acceleration.x, m_particle_acceleration.y, m_particle_acceleration.z);
        intersectParticlesGPU( TERRAIN_MIN_HEIGHT, halfDomain, m_dx, m_dxInv, m_VeffSplash, m_splashHeightChange, m_VeffSpray, m_VeffSplash );

        //copy back the changed data
        copybackGPU(DEPTH, m_depthField);
        copybackGPU(HEIGHT, m_heightField);

        //read back splash and spray particles
        copybackGPU(PARTICLE_POSITIONS, (float*) m_splash_positions);
        copybackGPU(PARTICLE_VELOCITIES, (float*) m_splash_velocities);
        copybackGPU(SPRAY_POSITIONS, (float*) m_spray_positions);
        copybackGPU(SPRAY_VELOCITIES, (float*) m_spray_velocities);
        copybackGPU(FOAM_POSITIONS, (float*) m_foam_positions);
        copybackGPU(FOAM_TTLS, (float*) m_foam_ttls);

        copybackGPU(VEL_U, m_velocityU);
        copybackGPU(VEL_W, m_velocityW);

        //read back foam information
        //add new foam particles
        copybackGPU(SPLASH_TO_FOAM, m_splash_to_foam_array);
        addNewFoamParticles();

        //check for breaking waves
        checkBreakingWavesGPU( ALPHA_MIN_SPLASH * GRAVITY * m_dt * m_dxInv, V_MIN_SPLASH, L_MIN_SPLASH,
                                    m_dx, m_dxInv, m_dt, halfDomain, m_splashHeightChange, LAMBDA_Y, BREAKING_WAVE_NUM_SPLASH_PARTICLES );
        copybackGPU( BREAKING_WAVES, m_breakingWavesGrid);
        addBreakingWaveParticles();
    }

    if(m_glw->m_useDampening) dampenWaves();

    m_timeElapsed += dt;
    m_updateCount++;

#ifdef SAVE_IMAGE
    const int savePerFrames = 5;
    if( (m_updateCount%savePerFrames) == 0 )
        saveToImage( VELOCITY );
#endif

}

void FluidGPU::addDrop(const int posX, const int posZ)
{
    // Fixed size
    int radius = m_gridSize/25;
    if( radius < 1 )
        radius = 1;
    float h = m_domainSize/15;

    addDropGPU( posX, posZ, radius, h );
    //copybackGPU(DEPTH,m_depthField);
    copybackGPU(HEIGHT,m_heightField);
    copybackGPU( PAINT, (float*)m_paintField );
}

/**
 * @brief Add drop to random positions
 *
 */
void FluidGPU::addRandomDrop( const float freq )
{
    /**
     * This function is not available
    **/
    // freq is not useful right now
    //float rnd = randomFloatGenerator();
    //int posX = rand()%m_gridSize;
    //int posY = rand()%m_gridSize;
    // int radius = rand()%(m_gridSize/15);
    /**
         * @brief not applicable right now
         */
    /**   int posX = 0.5*m_gridSize;
        int posZ = 0.5*m_gridSize;
        int radius = rand()%(m_gridSize/15);
    float rndHeight = randomFloatGenerator( 5, 15 );

    incrementH( posX, posZ, radius, rndHeight );**/
}

/**
 * @brief setColor Set the fluid's color
 * @param color The parameter
 */
void FluidGPU::setColor(const Colorf color)
{
    m_color = color;
}

void FluidGPU::backupHeight( Terrain* t )
{
    std::cout<<"Back up terrain heights into fluid..."<<std::endl;
    assert( t->getGridLength() == m_gridSize );

    // assert( t->getBound() == m_domainSize/2 );
    const Vector3* terrainVertices = t->getVerts();
    for( int i = 0; i < m_gridSize; i++ )
    {
        for( int j = 0; j < m_gridSize; j++ )
        {
            const int index = getIndex1D( i,j, TERRAINH);
            m_terrainHeightField[index] = terrainVertices[index].y;
        }
    }
    initBuffer();
    float halfDomain = m_domainSize/2.f;
    initGridGPU( m_gridSize, m_gridPaintSize, m_dx, halfDomain, m_terrainHeightField );
    copybackGPU(PAINT,(float*)m_paintField);
  //  copybackGPU(DEPTH,m_depthField);
    copybackGPU(HEIGHT,m_heightField);
    copybackGPU(NORMAL, (float*)m_paintNormalField );
    std::cout<<"Finshed backup"<<std::endl;

    //initialize dampening maps
    initDampeningFieldsGPU( DAMPENING_REGION, QUADRATIC_A, QUADRATIC_B, QUADRATIC_C );
    copybackGPU(SIGMA, m_sigmaField);
    copybackGPU(GAMMA, m_gammaField);
    copybackGPU(PHI, m_psiField);
    copybackGPU(PSI, m_phiField);

    //initialize particles
    initParticlesGPU(TERRAIN_MIN_HEIGHT, TOTAL_NUM_SPLASH_PARTICLES, TOTAL_NUM_SPRAY_PARTICLES, TOTAL_NUM_FOAM_PARTICLES);

    if(m_glw->m_useParticles)
    {
        copybackGPU(PARTICLE_POSITIONS, (float*)m_splash_positions );
        copybackGPU(PARTICLE_VELOCITIES, (float*)m_splash_velocities );
        copybackGPU(SPRAY_POSITIONS, (float*) m_spray_positions);
        copybackGPU(SPRAY_VELOCITIES, (float*) m_spray_velocities);
        copybackGPU(FOAM_POSITIONS, (float*) m_foam_positions);
        copybackGPU(FOAM_TTLS, (float*) m_foam_ttls);
    }

}

/**
 * @brief init Initialize the variables
 * @param gridSize The length of the grid
 */
void FluidGPU::init(const int gridSize, const float domainSize)
{
    m_gridSize = gridSize;
#ifdef RENDER_VOLUME
    m_gridPaintSize = gridSize + 2;
#else
    m_gridPaintSize = gridSize;
#endif

    m_uWidth = m_gridSize + 1;
    m_domainSize = domainSize;
    m_dx = m_domainSize/(float)m_gridSize;
    m_dxInv = 1.f/m_dx;

    /**
     * The height array should be (m_girdSize)x(m_gridSize)
     * The u array should be (m_gridSize)x(m_gridSize+1)
     * The v array should be (m_gridSize+1)x(m_gridSize)
     */

    /**
     *  Initialize 1D array
     **/
    int size = m_gridSize*m_gridSize*sizeof(float);
    m_depthField = (float*)malloc(size);
    m_terrainHeightField = (float*)malloc(size);
    m_sigmaField = (float*)malloc(size);
    m_gammaField = (float*)malloc(size);
    m_phiField = (float*)malloc(size);
    m_psiField = (float*)malloc(size);
    m_heightField = (float*)malloc(size);
    m_breakingWavesGrid = (float*)malloc(size);

    int paintSize = (m_gridPaintSize)*(m_gridPaintSize)*sizeof(Vector3);
    m_paintField = (Vector3*)malloc(paintSize);

    for( int i = 0; i < m_gridSize*m_gridSize; i++ )
    {
        m_depthField[i] = 0.f;
        m_terrainHeightField[i] = 0.f;
        m_sigmaField[i] = 0.f;
        m_gammaField[i] = 0.f;
        m_phiField[i] = 0.f;
        m_psiField[i] = 0.f;
        m_heightField[i] = 0.f;
        m_breakingWavesGrid[i] = 0.f;
    }

    size = (m_gridSize)*(m_gridSize+1)*sizeof(float);
    m_velocityU = (float*)malloc(size);
    m_velocityW = (float*)malloc(size);
    for( int i = 0; i < (m_gridSize)*(m_gridSize+1); i++ )
    {
        m_velocityU[i] = 0.f;
        m_velocityW[i] = 0.f;
    }

    m_paintNormalField = (Vector3*)malloc(paintSize);
    for( int i = 0; i < m_gridPaintSize*m_gridPaintSize; i++ )
    {
        // Initial value for normals
        m_paintNormalField[i] = Vector3(0,1,0);
    }

    int sizeIndex = (m_gridSize+1)*(m_gridSize-1)*2;
    m_indices = (GLuint*)malloc(sizeIndex*sizeof(sizeIndex));
    initIndices();
    buildTriangleList();
    // Set the random seed
    srand((unsigned)time(0));
    m_updateCount = 0;
    m_timeElapsed = 0.f;

    // Default color
    m_color = Colorf(0.1f,0.4f,0.8f,1.0f);
    m_renderNormals = false;

    // Find if cuda is supported
    if( !findSupportDevice() )
    {
        printf("GPU does not support CUDA!\n");
    }
#ifdef USE_PARTICLES
    initParticleSources();
#endif

    if(m_glw->m_useParticles)
    {
        //universal constants
        m_particle_acceleration = Vector3(0, GRAVITY, 0);

        //spray array
        int spraySize = sizeof(Vector3) * TOTAL_NUM_SPRAY_PARTICLES;
        m_spray_positions = (Vector3*) malloc(spraySize);
        m_spray_velocities = (Vector3*) malloc(spraySize);
        for(int i = 0; i < TOTAL_NUM_SPRAY_PARTICLES; i++){
            m_spray_positions[i] = Vector3(0, -1, 0);
            m_spray_velocities[i] = Vector3::zero();
        }
        m_last_active_spray_index = -1;

        //splash array
        int splashSize = sizeof(Vector3) * TOTAL_NUM_SPLASH_PARTICLES;
        m_splash_positions = (Vector3*) malloc(splashSize);
        m_splash_velocities = (Vector3*) malloc(splashSize);
        for(int i = 0; i < TOTAL_NUM_SPLASH_PARTICLES; i++){
            m_splash_positions[i] = Vector3(0, -1, 0);
            m_splash_velocities[i] = Vector3::zero();
        }
        m_last_active_splash_index = -1;

        //foam array
        int foamSize = sizeof(Vector3) * TOTAL_NUM_FOAM_PARTICLES;
        m_foam_positions = (Vector3*) malloc(foamSize);
        m_foam_ttls = (float*) malloc(foamSize);
        for(int i = 0; i < TOTAL_NUM_FOAM_PARTICLES; i++){
            m_foam_positions[i] = Vector3(0, -1, 0);
            m_foam_ttls[i] = 0.0f;
        }
        m_last_active_foam_index = -1;

        //splash to foam array
        m_splash_to_foam_array = (float*) malloc(splashSize);
        for(int i = 0; i < TOTAL_NUM_SPLASH_PARTICLES; i++){
            m_splash_to_foam_array[i] = -1;
        }

        //pre-caclulating several needed values to speed up things throughout
        m_VeffSpray = C_DEPOSIT * (4.0 / 3.0) * M_PI *
                SPRAY_PARTICLE_RADIUS * SPRAY_PARTICLE_RADIUS * SPRAY_PARTICLE_RADIUS;
        m_sprayHeightChange = m_VeffSpray / (m_dx * m_dx);

        m_VeffSplash = C_DEPOSIT * (4.0 / 3.0) * M_PI *
                SPLASH_PARTICLE_RADIUS * SPLASH_PARTICLE_RADIUS * SPLASH_PARTICLE_RADIUS;
        m_splashHeightChange = m_VeffSplash / (m_dx * m_dx);


    }
}

/**
 * @brief initBuffer Initialize the buffer
 */
void FluidGPU::initBuffer()
{
    assert( m_paintField != NULL );
    int terrainSize;

    terrainSize = (m_gridPaintSize)*(m_gridPaintSize);

    glGenBuffers(1,&m_vertexBuffer );
    glBindBuffer( GL_ARRAY_BUFFER, m_vertexBuffer );
    glBufferData( GL_ARRAY_BUFFER, terrainSize*sizeof(Vector3), m_paintNormalField, GL_STATIC_DRAW );

    glGenBuffers( 1, &m_normalBuffer );
    glBindBuffer( GL_ARRAY_BUFFER, m_normalBuffer );
    glBufferData( GL_ARRAY_BUFFER, terrainSize*sizeof(Vector3), m_paintNormalField, GL_STATIC_DRAW );

    int indexSize;
    indexSize = (m_gridPaintSize+1)*(m_gridPaintSize-1)*2;

    glGenBuffers(1,&m_indexBuffer );
    glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_indexBuffer );
    glBufferData( GL_ELEMENT_ARRAY_BUFFER, indexSize*sizeof(GLuint), m_indices, GL_STATIC_DRAW );
}

/**
 * @brief init index array
 */
void FluidGPU::initIndices()
{
    assert(m_indices != NULL );
    int index = 0;
    for( int z = 0; z < m_gridPaintSize-1; z++ )
    {
        int x;
        for( x =  0; x < m_gridPaintSize; x++ )
        {

            m_indices[index] = x + z*m_gridPaintSize;
            index++;
            m_indices[index] = x + (z+1)*m_gridPaintSize;
            index++;
        }
        m_indices[index] = x-1 + (z+1)*m_gridPaintSize;
        index++;
        m_indices[index] = 0 + (z+1)*m_gridPaintSize;
        index++;
    }
}

/**
 * @brief Check boundary
 */
void FluidGPU::checkBoundary()
{
}

/**
 * @brief Write the height field or velocity to image
 */
void FluidGPU::saveToImage( FieldType type )
{
    QImage saveImage = QImage(m_gridSize,m_gridSize,QImage::Format_RGB32);
    BGRA* data =(BGRA*)saveImage.bits();
    float maxH = maxHeight*1.2;
    float minH = defaultHeight*0.5;
    float minV = -2;
    float maxV = 2;
    for( int i = 0; i < m_gridSize; i++ )
    {
        for( int j = 0; j < m_gridSize; j++ )
        {
            int index = getIndex1D(i,j,DEPTH);
            float curH = m_depthField[index];
            int gray = (int)((curH - minH)/(maxH - minH)*255);
            data[index].r = gray;
            data[index].g = gray;
            data[index].b = gray;
            data[index].a = 255;
            if( type == VEL )
            {
                float u = 0.5*(m_velocityU[getIndex1D(i,j,VEL_U)] + m_velocityU[getIndex1D(i,j+1,VEL_U)]);
                float w = 0.5*(m_velocityW[getIndex1D(i,j,VEL_W)] + m_velocityW[getIndex1D(i+1,j,VEL_W)]);
                int iu = (int)(((u - minV)/(maxV - minV))*100.f);
                int iw = (int)(((w - minV)/(maxV - minV))*100.f);
                data[index].r +=  iu;
                data[index].g += iw;
                if( data[index].r > 255 )
                    data[index].r = 255;
                if( data[index].g > 255 )
                    data[index].b = 255;
            }
        }
    }
    std::stringstream ss;
    std::string fileName;
    std::string form = ".png";
    if( type == DEPTH )
    {
        std::string name = SAVE_NAME_DEPTH;
        ss<<name<<m_updateCount<<form;
        fileName = ss.str();
    }
    else if( type == VEL )
    {
        std::string name = SAVE_NAME_VELOCITY;
        ss<<name<<m_updateCount<<form;
        fileName = ss.str();
    }
    else
    {
        assert(0);
    }
    assert( fileName.size() > 0 );
    std::cout << "Save picture into "<<fileName<<std::endl;
    saveImage.save( fileName.c_str() );
}

/**
 * @brief drawFluid Draw the fluid with different method
 * @param method The method for drawing, could be DRAW_POINTS or DRAW_MESH
 */
void FluidGPU::drawFluid( DrawMethod method ) const
{
    if( method == DRAW_POINTS )
    {
        // Not applicable
        /*glPushMatrix();
        glBegin(GL_POINTS);
        glColor4f( m_color.r,m_color.g, m_color.b, m_color.a );
        for( int i = 0; i < m_gridPaintSize; i++ )
        {
            for( int j =0; j < m_gridPaintSize; j++ )
            {
                glVertex3fv(m_paintField[getIndex1D(i,j,PAINT)].xyz);
            }
        }
        glPopMatrix();
        glEnd();*/
    }
    else if ( method == DRAW_MESH_STRIP )
    {
        // This mode doesn't deal with the invisible triangles
        for( int i = 0; i < m_gridPaintSize-1; i++ )
        {
            glBegin( GL_TRIANGLE_STRIP );
            glColor4f(m_color.r,m_color.g,m_color.b,m_color.a);
            for( int j = 0; j < m_gridPaintSize; j++ )
            {
                int index = i*m_gridPaintSize + j;
                glNormal3fv( m_paintNormalField[index].xyz);
                glVertex3fv( m_paintField[index].xyz);
                index = (i+1)*m_gridPaintSize + j;
                glNormal3fv( m_paintNormalField[index].xyz );
                glVertex3fv( m_paintField[index].xyz);
            }
            glEnd();
        }
    }
    else if( method == DRAW_MESH )
    {
        // Currently disabled
//        float halfDomain = m_domainSize/2;
//        // Else we use triangles mode to draw(This will hide the invisible triangles)
//        glBegin( GL_TRIANGLES );
//        glColor4f(m_color.r,m_color.g,m_color.b,m_color.a);
//        for( int i = 0;  i < m_triangles.size(); i++ )
//        {
//            Tri curTri  = m_triangles[i];
//            int r[3]; int c[3];
//            r[0] = curTri.a2D.indRow; c[0] = curTri.a2D.indCol;
//            r[1] = curTri.b2D.indRow; c[1] = curTri.b2D.indCol;
//            r[2] = curTri.c2D.indRow; c[2] = curTri.c2D.indCol;
//            // Firstly check if the triangle is visible
//            /*              int count = 0;

//                for( int m = 0; m < 3; m++ )
//                {
//                    if( m_depthField[getIndex1D(r[m],c[m],DEPTH)] > EPSILON )
//                        count++;
//                }
//                if( count == 0 )
//                    continue;
//*/
//            float tx, ty, tz;
//            for( int m = 0; m < 3; m++ )
//            {
//                //if( m_depthField[getIndex1D(r[m],c[m],DEPTH)] > EPSILON )
//                glColor4f(m_color.r,m_color.g,m_color.b,m_color.a);
//                /*else
//                        glColor4f(0.f,0.f,0.f, m_color.a );
//                    */
//                int index = r[m]*m_gridPaintSize + c[m];
//                tx = -halfDomain + c[m]*m_dx;
//                ty = m_heightField[index];
//                if( m_depthField[index] < 0.5 )
//                    ty = ty - 0.5;
//                tz = - halfDomain + r[m]*m_dx;
//                glNormal3f( m_paintNormalField[index].x, m_paintNormalField[index].y, m_paintNormalField[index].z );
//                glVertex3f( tx, ty, tz );
//            }
//        }
//        glEnd();
    }
    else if( method = DRAW_MESH_VBO )
    {
        glBindBuffer( GL_ARRAY_BUFFER, m_vertexBuffer );
        Vector3* vertBuffer = (Vector3*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY );

        memcpy( vertBuffer, m_paintField, sizeof(Vector3)*(m_gridPaintSize)*(m_gridPaintSize));

        glUnmapBuffer(GL_ARRAY_BUFFER);

         glBindBuffer( GL_ARRAY_BUFFER, m_vertexBuffer );
        glVertexPointer(3,GL_FLOAT,0,(char*)NULL);
        glEnableClientState( GL_VERTEX_ARRAY );

        glBindBuffer( GL_ARRAY_BUFFER, m_normalBuffer );
        Vector3* normBuffer = (Vector3*)glMapBuffer(GL_ARRAY_BUFFER, GL_WRITE_ONLY );
        memcpy( normBuffer, m_paintNormalField, sizeof(Vector3)*(m_gridPaintSize)*(m_gridPaintSize));
        glUnmapBuffer(GL_ARRAY_BUFFER);

        glBindBuffer( GL_ARRAY_BUFFER, m_normalBuffer );
        glNormalPointer(GL_FLOAT,0,(char*)NULL);
        glEnableClientState( GL_NORMAL_ARRAY );

        glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, m_indexBuffer );

        int indexSize;
        indexSize = (m_gridPaintSize + 1)*(m_gridPaintSize- 1)*2;
        glColor4f(m_color.r,m_color.g,m_color.b,m_color.a);

        glDrawElements( GL_TRIANGLE_STRIP, indexSize, GL_UNSIGNED_INT, 0 );

        glDisableClientState( GL_NORMAL_ARRAY );
        glDisableClientState( GL_VERTEX_ARRAY );
        glBindBuffer( GL_ARRAY_BUFFER, 0 );
         glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, 0 );
    }
    else
    {
        assert(0);
    }
}

/**
 * @brief Draw the normals of the fluid points
 */
void FluidGPU::drawNormal() const
{
    if (m_renderNormals)
    {
        glColor3f(1,1,1);
        const int magn = 2;
        for (int row = 0; row < m_gridPaintSize; row++)
        {
            for (int column = 0; column <  m_gridPaintSize; column++)
            {
                glBegin(GL_LINES);

                Vector3 curVert = m_paintField[getIndex1D(row,column,PAINT)];
                Vector3 curNorm = m_paintNormalField[getIndex1D(row,column,NORMAL)];

                glNormal3f(curNorm.x,curNorm.y,curNorm.z);
                glVertex3f(curVert.x, curVert.y, curVert.z);
                glVertex3f(curVert.x +magn*curNorm.x,
                           curVert.y + magn*curNorm.y,
                           curVert.z + magn*curNorm.z);

                glEnd();
            }
        }
    }
}

/**
 * @brief build the triangle List
 */
void FluidGPU::buildTriangleList()
{
    cout<<"Build the triangle list..."<<endl;
    for( int i = 0; i < m_gridSize-1; i++ )
    {
        for( int j = 0; j < m_gridSize-1; j++ )
        {
            TriIndex t1,t2;
            t1.a2D.indRow = i;
            t1.a2D.indCol = j;
            t1.b2D.indRow = i+1;
            t1.b2D.indCol = j;
            t1.c2D.indRow = i;
            t1.c2D.indCol = j+1;
            m_triangles.push_back(t1);

            t2.a2D.indRow = i;
            t2.a2D.indCol = j+1;
            t2.b2D.indRow = i+1;
            t2.b2D.indCol = j;
            t2.c2D.indRow = i+1;
            t2.c2D.indCol = j+1;
            m_triangles.push_back(t2);
        }
    }
    cout<<"Finish building the triangle list."<<endl;
}

/**
 * @brief dampen waves for open water scenes
 */
void FluidGPU::dampenWaves(){
    //compute hRest
    float hRest = computeHRest();

    //dampen waves on GPU
    dampenWavesGPU( hRest, m_dt, m_dxInv, LAMBDA_UPDATE, LAMBDA_DECAY );

    //copyback
    copybackGPU(DEPTH, m_depthField);
    copybackGPU(HEIGHT, m_heightField);
    //copybackGPU(VEL_U, m_velocityU);
    //copybackGPU(VEL_W, m_velocityU);

    //TODO: needed?
    //copybackGPU(SIGMA, m_sigmaField);
    //copybackGPU(GAMMA, m_gammaField);
    //copybackGPU(PHI, m_phiField);
    //copybackGPU(PSI, m_psiField);
}

/**
 * @brief computes the resting height of the fluid
 * the resting height is computed as the average across the depth field
 * the result of this function is used in dampenWaves()
 * @return h_rest
 */
float FluidGPU::computeHRest(){
    float hRest = 0;
    float denominator = 0;
    for(int i = 0; i < m_gridSize * m_gridSize; i++){
        if(m_depthField[i] > 0){
            hRest += m_heightField[i];
            denominator += 1;
        }
    }
    return (hRest / denominator);
}

void FluidGPU::clampFields(){
    // clamp h(i,j) >= 0
    // clamp u(i,j) < alpha * (delta_x / delta_t)
    // clamp w(i,j) < alpha * (delta_x / delta_t)

    float velocityClamp = CLAMP_ALPHA * (m_dx / m_dt);
    clampFieldsGPU( velocityClamp );
}

/**
 * @brief getIndex1D Return the corresponding 1D index based on which type
 * @param i The row number
 * @param j The col number
 * @param type The type of the field
 * @return The 1D index
 */
int FluidGPU::getIndex1D( int i, int j, FieldType type) const
{
    switch( type )
    {
    case PAINT:
    case NORMAL:
        return i*m_gridPaintSize + j;
        break;
    case HEIGHT:
    case TERRAINH:
    case VEL_W:
    case DEPTH:
        return i*m_gridSize + j;
        break;
    case TMP:
    case VEL_U:
    {
        return i*(m_uWidth)+j;
        break;
    }
    default:
        printf("Invalid type\n");
        return -1;
    }
}

/**
 * @brief getFieldArray Get the pointer to the array with specified type
 * @param buffLength The buffer length returned
 * @return The buffer
 */
float* FluidGPU::getFieldArray( FieldType type, int& buffLength ) const
{
    switch( type )
    {
    case HEIGHT:
        buffLength = m_gridSize*m_gridSize;
        return m_heightField;
        break;
    case TERRAINH:
        buffLength = m_gridSize*m_gridSize;
        return m_terrainHeightField;
        break;
    case VEL_W:
        buffLength = m_gridSize*(m_gridSize+1);
        return m_velocityW;
        break;
    case DEPTH:
        buffLength = m_gridSize*m_gridSize;
        return m_depthField;
        break;
    case VEL_U:
    {
        buffLength = (m_gridSize+1)*m_gridSize;
        return m_velocityU;
        break;
    }
    default:
        printf("Invalid type\n");
        return 0;
    }
}

void FluidGPU::updateParticles(){
    //advance current particles
    for(int i = 0; i < m_particles.size(); i++){
        Particle *currParticle = m_particles[i];
        currParticle->updateParticle(m_dt);
    }

    //remove particles
    removeParticles();

    //add any new particles
    //checkForBreakingWaves();

    //copy back changes
   // copybackGPU(DEPTH,m_depthField);
   // copybackGPU(HEIGHT,m_heightField);
    copybackGPU(PAINT,(float*)m_paintField);

}

void FluidGPU::removeParticles(){
    double halfDomain = m_domainSize / 2.0;

    QVector<Vector3> values;

    int index = 0;
    while(index < m_particles.size()){
        //get current particle
        Particle *currParticle = m_particles[index];
        Vector3 position = currParticle->getPosition();

        //check bounds
        //check return to fluid
        if(position.x < -halfDomain || position.x > halfDomain ||
                position.z < -halfDomain || position.z > halfDomain ||
                position.y < 0){
            m_particles.remove(index);
            delete currParticle;
        } else {
            Vector3 currValues = fluidParticleInteractionCheck(currParticle);

            if(currValues.x >= 0){
                bool added = false;
                for(int i = 0; i < values.size(); i++){
                    if(currValues.x == values[i].x && currValues.z == values[i].z){
                        values[i].y += currValues.y;
                        added = true;
                        break;
                    }
                }
                if(!added){
                    values.append(currValues);
                }

                m_particles.remove(index);
                delete currParticle;
            } else {
                index++;
            }
        }
    }

    if(values.size() > 0){
        fluidParticleUpdate(values);
    }
}

bool FluidGPU::fluidParticleInteraction(Particle *particle){
    Vector3 position = particle->getPosition();

    //get positions
    float halfDomain = m_domainSize / 2.0;
    float lenX = position.x + halfDomain;
    float lenZ = position.z + halfDomain;

    int j = (int) min(m_gridSize, max(0.0, round(lenX / m_dx)));
    int i = (int) min(m_gridSize, max(0.0, round(lenZ / m_dx)));

    int depthIndex = getIndex1D(i, j, DEPTH);
    int heightIndex = getIndex1D(i, j, HEIGHT);
    int velocityUIndex = getIndex1D(i, j, VEL_U);
    int velocityWIndex = getIndex1D(i, j, VEL_W);

    if(position.y <= m_heightField[heightIndex]){
        float Veff = particle->getVolume();
        float heightChange = Veff / (m_dx * m_dx);

        //add to height
        addDropGPU( j, i, 2, heightChange );

        float denominator = (m_depthField[depthIndex] * m_dx * m_dx) + Veff;

        //TODO: make this do something
        m_velocityU[velocityUIndex] = ((m_velocityU[depthIndex] * m_depthField[depthIndex] * m_dx * m_dx) + (particle->getVelocity().x * Veff)) / denominator;
        m_velocityW[velocityWIndex] = ((m_velocityW[depthIndex] * m_depthField[depthIndex] * m_dx * m_dx) + (particle->getVelocity().z * Veff)) / denominator;

        return true;
    }
    return false;
}

Vector3 FluidGPU::fluidParticleInteractionCheck(Particle *particle){
    Vector3 position = particle->getPosition();

    //get positions
    float halfDomain = m_domainSize / 2.0;
    float lenX = (position.x + halfDomain) * m_dxInv;
    float lenZ = (position.z + halfDomain) * m_dxInv;

    int j = (int) min(m_gridSize, max(0.0, round(lenX)));
    int i = (int) min(m_gridSize, max(0.0, round(lenZ)));

    int heightIndex = getIndex1D(i, j, HEIGHT);

    if(position.y <= m_heightField[heightIndex]){
        float heightChange = particle->getVolume() / (m_dx * m_dx);
        return Vector3((float)j, heightChange, (float)i);
    }
    return Vector3(-1, 0, -1);
}

void FluidGPU::fluidParticleUpdate(QVector<Vector3> values){
    for(int p = 0; p < values.size(); p++){
        int j = (int)values[p].x;
        int i = (int)values[p].z;
        float heightChange = values[p].y;
        addDropGPU( j, i, 2, heightChange );
    }
}

void FluidGPU::drawParticles() const{
    //begin
    //glBegin(GL_QUADS);
    glEnable(GL_PROGRAM_POINT_SIZE_EXT);
    glPointSize(20);
    glBegin(GL_POINTS);
    glEnable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_ONE,GL_ONE);

    //iterate through and draw particles as points
    for(int i = 0; i < m_particles.size(); i++){
        Particle *drawParticle = m_particles[i];
        drawParticle->drawParticle();
    }

    //end
    glDisable(GL_BLEND);
    glDisable(GL_CULL_FACE);
    //glDisable(GL_PROGRAM_POINT_SIZE_EXT);
    glEnd();
}

void FluidGPU::addDroppingParticles(const int posX, const int posZ){
    //support values
    float radius = m_dx * PARTICLE_DROPPING_RADIUS;
    float Veff = C_DEPOSIT * (4 / 3) * M_PI *
            SPLASH_PARTICLE_RADIUS * SPLASH_PARTICLE_RADIUS * SPLASH_PARTICLE_RADIUS;
    float halfDomain = m_domainSize / 2.0;

    float coordX = -halfDomain + (posX * m_dx);
    float coordY = -halfDomain + (posZ * m_dx);

#ifdef USE_PARTICLES
//    for(int i = 0; i < NUM_DROPPING_PARTICLES; i++){
//        //jitter particle positions
//        float randX = randomFloatGenerator(-radius, radius);
//        float randY = randomFloatGenerator(0, PARTICLE_DROP_RANGE);
//        float randZ = randomFloatGenerator(-radius, radius);

//        Vector3 position = Vector3(coordX + randX, PARTICLE_DROP_HEIGHT + randY, coordY + randZ);
//        Vector3 velocity = Vector3::zero();
//        Vector3 acceleration = Vector3(0, GRAVITY, 0);

//        //make new particle
//        Particle *newParticle = new Particle(SPLASH_PARTICLE_RADIUS, Veff, position, velocity, acceleration);
//        m_particles.append(newParticle);
//    }
#endif

    if(m_glw->m_useParticles)
    {
        int currIndex = 0;
        for(int i = 0; i < NUM_DROPPING_PARTICLES; i++){
            if(currIndex >= TOTAL_NUM_SPRAY_PARTICLES){
                break;
            }

            //jitter particle positions
            float randX = randomFloatGenerator(-radius, radius);
            float randY = randomFloatGenerator(0, PARTICLE_DROP_RANGE);
            float randZ = randomFloatGenerator(-radius, radius);

            Vector3 position = Vector3(coordX + randX, PARTICLE_DROP_HEIGHT + randY, coordY + randZ);
            Vector3 velocity = Vector3::zero();

            //find new inactive particle
            bool added = false;
            while(currIndex < TOTAL_NUM_SPRAY_PARTICLES && !added){
                //if(m_spray_positions[currIndex].y < 0){
                if(m_spray_positions[currIndex].y < TERRAIN_MIN_HEIGHT){
                    //set new particle as active
                    m_spray_positions[currIndex] = position;
                    m_spray_velocities[currIndex] = velocity;

                    added = true;
                    break;
                }

                currIndex++;
            }
        }

        inputParticlesGPU((float*)m_spray_positions, (float*)m_spray_velocities);
    }

}

void FluidGPU::initParticleSources(){
    float halfDomain = m_domainSize / 2.0;

    //    //make a big, low flow source
    //    Vector3 startingCorner = Vector3(-halfDomain, 60, -halfDomain);
    //    Vector3 endingCorner = Vector3(halfDomain, 65, halfDomain);
    //    ParticleSource *particleSource = new ParticleSource(startingCorner, endingCorner, 50);
    //    m_particle_sources.append(particleSource);

    //    //make a small, high flow source
    //    Vector3 startingCorner2 = Vector3(-m_dx, 100, -m_dx);
    //    Vector3 endingCorner2 = Vector3(m_dx, 110, m_dx);
    //    Vector3 startingCorner2 = Vector3(-halfDomain, 100, -5);
    //    Vector3 endingCorner2 = Vector3(halfDomain, 120, 5);
    //    ParticleSource *particleSource2 = new ParticleSource(startingCorner2, endingCorner2, 100);
    //    m_particle_sources.append(particleSource2);
}

void FluidGPU::updateParticleSources(){
    for(int i = 0; i < m_particle_sources.size(); i++){
        QVector<Particle*> newParticles = m_particle_sources[i]->generateParticles();

        for(int j = 0; j < newParticles.size(); j++){
            m_particles.append(newParticles[j]);
        }
    }
}

void FluidGPU::drawParticles2() {
    if(m_glw->m_useParticles)
    {
        //begin
        //drawings spray only
        //@NOTE these next two lines must come outside the glBegin/glEnd block
//        glEnable(GL_PROGRAM_POINT_SIZE_EXT);
//        glPointSize(0.5);
        glColor4f(SPRAY_COLOR);
        glBegin(GL_POINTS);
//        glEnable(GL_CULL_FACE);
//        glEnable(GL_BLEND);
//        glBlendFunc(GL_ONE,GL_ONE);
        //iterate through each position and draw a point
        for(int i = 0; i < TOTAL_NUM_SPRAY_PARTICLES; i++){
            Vector3 position = m_spray_positions[i];
            //if(position.y >= 0){
            if(position.y >= TERRAIN_MIN_HEIGHT){
                glVertex3f(position.x, position.y, position.z);
            }
        }

        //end
//        glDisable(GL_BLEND);
//        glDisable(GL_CULL_FACE);
        glEnd();

        //now draw splash a different size and color
//        glEnable(GL_PROGRAM_POINT_SIZE_EXT);
//        glPointSize(1.5);
        glColor4f(SPLASH_COLOR);
        glBegin(GL_POINTS);
//        glEnable(GL_CULL_FACE);
//        glEnable(GL_BLEND);
//        glBlendFunc(GL_ONE,GL_ONE);

        //iterate through each position and draw a point
        for(int i = 0; i < TOTAL_NUM_SPLASH_PARTICLES; i++){
            Vector3 position = m_splash_positions[i];
            //if(position.y >= 0){
            if(position.y >= TERRAIN_MIN_HEIGHT){
                glVertex3f(position.x, position.y, position.z);
            }
        }

        //end
//        glDisable(GL_BLEND);
//        glDisable(GL_CULL_FACE);
        glEnd();


        //now draw foam a different size and color
//        glEnable(GL_PROGRAM_POINT_SIZE_EXT);
//        glPointSize(2);
        glColor4f(FOAM_COLOR);
        glBegin(GL_POINTS);
//        glEnable(GL_CULL_FACE);
//        glEnable(GL_BLEND);
//        glBlendFunc(GL_ONE,GL_ONE);

        //iterate through each position and draw a point
        for(int i = 0; i < TOTAL_NUM_FOAM_PARTICLES; i++){
            Vector3 position = m_foam_positions[i];
            //if(position.y >= 0){
            if(position.y >= TERRAIN_MIN_HEIGHT){
                glVertex3f(position.x, position.y, position.z);
            }
        }

        //end
        //        glDisable(GL_BLEND);
        //        glDisable(GL_CULL_FACE);
        glEnd();


    }
}

void FluidGPU::addBreakingWaveParticles(){
    if(m_glw->m_useParticles)
    {
        float halfDomain = m_domainSize / 2.0;
        float halfDx = m_dx / 2.0;

        int currSplashIndex = 0;
        int currSprayIndex = 0;
        for(int i = 1; i < m_gridSize - 1; i++){
            if(currSplashIndex >= TOTAL_NUM_SPLASH_PARTICLES && currSprayIndex >= TOTAL_NUM_SPRAY_PARTICLES){
                break;
            }

            for(int j = 1; j < m_gridSize - 1; j++){
                if(currSplashIndex >= TOTAL_NUM_SPLASH_PARTICLES && currSprayIndex >= TOTAL_NUM_SPRAY_PARTICLES){
                    break;
                }

                //index of current cell
                int index = getIndex1D(i, j, DEPTH);
                int uindex = getIndex1D(i, j, VEL_U);
                int windex = getIndex1D(i, j, VEL_W);

                int numToAdd = (int) floor(m_breakingWavesGrid[index]);

                //TODO: position of current grid cell
                float posX = -halfDomain + (j * m_dx);
                float posZ = -halfDomain + (i * m_dx);

                for(int a = 0; a < numToAdd; a++){
                    //jitter particle positions
                    float randX = randomFloatGenerator(-halfDx, halfDx);
                    float randZ = randomFloatGenerator(-halfDx, halfDx);

                    Vector3 position = Vector3(posX + randX, m_heightField[index], posZ + randZ);
                    //TODO: wrong y component!
                    Vector3 velocity = Vector3(BREAKING_WAVE_VEL_MULTIPLIER * m_velocityU[uindex],
                                               BREAKING_WAVE_VEL_MULTIPLIER * LAMBDA_Y * V_MIN_SPLASH,
                                               BREAKING_WAVE_VEL_MULTIPLIER * m_velocityW[windex]);

                    //if possible, add new splash particle
                    if(m_depthField[index] >= m_splashHeightChange){
                        //add new active splash particle
                        bool added = false;
                        while(currSplashIndex < TOTAL_NUM_SPLASH_PARTICLES && !added){
                            //if(m_splash_positions[currIndex].y < 0){
                            if(m_splash_positions[currSplashIndex].y < TERRAIN_MIN_HEIGHT){
                                //set new particle as active
                                m_splash_positions[currSplashIndex] = position;
                                m_splash_velocities[currSplashIndex] = velocity;

                                m_depthField[index] -= m_splashHeightChange;

                                added = true;
                                break;
                            }

                            currSplashIndex++;
                        }
                    }

                    //if possible, add new spray particle
                    if(m_depthField[index] >= m_sprayHeightChange){
                        //add new active spray particle
                        bool added = false;
                        while(currSprayIndex < TOTAL_NUM_SPRAY_PARTICLES && !added){
                            if(m_spray_positions[currSprayIndex].y < TERRAIN_MIN_HEIGHT){
                                //set new particle as active
                                m_spray_positions[currSprayIndex] = position;
                                m_spray_velocities[currSprayIndex] = SPRAY_VEL_MULTIPLIER * velocity;

                                m_depthField[index] -= m_sprayHeightChange;

                                added = true;
                                break;
                            }

                            currSprayIndex++;
                        }
                    }
                }
            }
        }

        //send particle information to CUDA
        inputParticlesGPU((float*)m_splash_positions, (float*)m_splash_velocities);
        inputSprayParticlesGPU((float*)m_spray_positions, (float*)m_spray_velocities);
        inputDepthGPU( m_depthField );
        copybackGPU(HEIGHT, m_heightField );
    }
}

void FluidGPU::addNewFoamParticles(){
    if(m_glw->m_useParticles)
    {
        float halfDomain = m_domainSize / 2.0;

        int currFoamIndex = 0;

        //iterate through splash particles
        for(int p = 0; p < TOTAL_NUM_SPLASH_PARTICLES; p++){
            if(currFoamIndex >= TOTAL_NUM_FOAM_PARTICLES){
                break;
            }

            //if this has just turned to foam
            if(m_splash_to_foam_array[p] > 0){
                //determine where the particle is
                Vector3 currPosition = m_splash_positions[p];
                float lenX = (currPosition.x + halfDomain) * m_dxInv;
                float lenZ = (currPosition.z + halfDomain) * m_dxInv;
                int i = (int) min(m_gridSize - 1, max(0.0, round(lenX)));
                int j = (int) min(m_gridSize - 1, max(0.0, round(lenZ)));
                int gridIndex = getIndex1D(i, j, HEIGHT);

                //get the position
                Vector3 position = m_splash_positions[p];
                position.y = m_heightField[gridIndex];

                //find a foam particle to activate, if possible
                bool added = false;
                while(currFoamIndex < TOTAL_NUM_FOAM_PARTICLES && !added){
                    if(m_foam_positions[currFoamIndex].y < TERRAIN_MIN_HEIGHT){
                        //set new particle as active
                        m_foam_positions[currFoamIndex] = position;
                        m_foam_ttls[currFoamIndex] = FOAM_TTL;

                        added = true;
                        break;
                    }

                    currFoamIndex++;
                }
            }
        }

        //send particle information to CUDA
        inputFoamParticlesGPU((float*)m_foam_positions, m_foam_ttls);
    }
}
