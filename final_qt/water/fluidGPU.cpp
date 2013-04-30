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
void initGridGPU( const int hostGridSize, const float hostdx, const float* hostTerrainMap );
void copybackGPU(FieldType type, float* hostMap  );
void destroyGPUmem();
void addDropGPU(const int posX, const int posZ, const int radius, const float h );
void advectGPU(const float dt);
void updateFluidGPU( const float dt );
bool findSupportDevice();
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

FluidGPU::FluidGPU( Terrain *t )
{
    init( t->getGridLength(), 2*t->getBound() );
}

FluidGPU::~FluidGPU()
{
    // Release the heap
    safeFreeArray1D( m_velocityU );
    safeFreeArray1D( m_velocityW );
    safeFreeArray1D( m_terrainHeightField );
    safeFreeArray1D( m_sigmaField );
    safeFreeArray1D( m_gammaField );
    safeFreeArray1D( m_phiField );
    safeFreeArray1D( m_psiField );
    safeFreeArray1D( m_normalField );
    safeFreeArray1D( m_heightField );
    safeFreeArray1D( m_depthField );

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

//Main algo loop - once per time step

    //Height field fluid simulation - Sec 2.1
        //discretize the simulation domain where the heights are stored at the cell centers
        //and the velocity components on faces
    //glPolygonMode(GL_FRONT, GL_LINE);
    drawFluid( DRAW_MESH );
#ifdef USE_PARTICLES
    drawParticles();
#endif
  if( m_renderNormals )
        drawNormal();
//  glPolygonMode(GL_FRONT,GL_FILL);
        //employing time splitting by first solving for self advection of the velocity field
            //advection is the unconditionally stable modified MacCormack method
            //but we fall back onto the semi-Lagrangian method if the first equation result is out of bounds

        //then integrating the height field and velocity field forward in time
            //explicity intergrate the height of the field using equation (7) to guarentee mass preservation
            //@TODO - modify this step to take waterfall discontinuties into account for Sec 2.4.3

            //update face velocities taking the gradient of the water height into account
            //@TODO - this also needs to be modified to account for waterfalls, Sec 2.4.3

        //boundary conditions
            //dealing with reflective and static surfaces by setting face velocity to zero
            //the entire "wet" cell needs to be higher than the terrain level in order to flow, if not treat as stopped

            //for open water scenes implement Perfectly Matched Layers to dampen out the waves when the get near the edge
            //the width of the dampening region is 10 cells

            //@NOTE: for stability, be sure to clamp h_i,j to always be >= zero

            //@NOTE:  for violent wave stabilty clamp u_i +1/2, j and v_i, j+1/2 to alpha dx/dt = 0.5

            //@NOTE: limit the water depth used for the height intergration, beta = 2

        //Overshooting reduction - Sec 2.2
        //to prevent triangle overlap when moving to shallow regions from deeper water
        //to detect edges of the waves and reduce the magnitude - need to fix in x and z planes

    //Solids simulation
    //we currently have the terrain right in glwidget - will either need to move that here to pass a pointer
    //for the next step of the algo
    //Two-way coupling of height field and solids - Sec 2.3

        //recursively divide each triangle in fluid and solids into sub triangles smaller than kdx^2 where k = 1
        //let p = the position at this time step and v = the velocity at this time step of the centroid of the sub triangle
        //A is the area of the subtriangle with  p and v obtained by baycentric interpolation from the original triangle
        //with normal of n from the vector triangle

        //calculate the buoyancy of the triangle
        //calculate the lift of the triangle
        //calculate the drag of the triangle

        //add these forces to the solid the subtriangle belongs to, for a fluid weight the forces to the three vertices of the original tri


        //modify the height and velocity of the fluid due to solids using algo 2


    //Particles generation and simulation - Sec 2.4

    //spray/splash and rform
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
    copybackGPU(DEPTH,m_depthField);
    copybackGPU(HEIGHT,m_heightField);
    copybackGPU(NORMAL,(float*)m_normalField);
    //clampFields();

#ifdef USE_PARTICLES
    updateParticles();
    updateParticleSources();
#endif

#ifdef DAMPEN_WAVES
    //dampenWaves();
#endif

//    computeNormal();

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
    int radius = m_gridSize/20;
    if( radius < 1 )
        radius = 1;
    float h = 4;

    addDropGPU( posX, posZ, radius, h );
    copybackGPU(DEPTH,m_depthField);
    copybackGPU(HEIGHT,m_heightField);
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
   initGridGPU( m_gridSize, m_dx, m_terrainHeightField );
   copybackGPU(DEPTH,m_depthField);
   copybackGPU(HEIGHT,m_heightField);
   copybackGPU(NORMAL, (float*)m_normalField );
   std::cout<<"Finshed backup"<<std::endl;
}

/**
 * @brief init Initialize the variables
 * @param gridSize The length of the grid
 */
void FluidGPU::init(const int gridSize, const float domainSize)
{
    m_gridSize = gridSize;
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
    for( int i = 0; i < m_gridSize*m_gridSize; i++ )
    {
        m_depthField[i] = 0.f;
        m_terrainHeightField[i] = 0.f;
        m_sigmaField[i] = 0.f;
        m_gammaField[i] = 0.f;
        m_phiField[i] = 0.f;
        m_psiField[i] = 0.f;
        m_heightField[i] = 0.f;
    }

    size = (m_gridSize)*(m_gridSize+1)*sizeof(float);
    m_velocityU = (float*)malloc(size);
    m_velocityW = (float*)malloc(size);
    for( int i = 0; i < (m_gridSize)*(m_gridSize+1); i++ )
    {
        m_velocityU[i] = 0.f;
        m_velocityW[i] = 0.f;
    }
    size = m_gridSize*m_gridSize*sizeof(Vector3);
    m_normalField = (Vector3*)malloc(size);
    for( int i = 0; i < m_gridSize*m_gridSize; i++ )
    {
        // Initial value for normals
        m_normalField[i] = Vector3(1,1,0);
    }
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

    // initializing wave dampening fields
    /*
    m_sigmaField.resize(m_gridSize);
    m_gammaField.resize(m_gridSize + 1);
    m_phiField.resize(m_gridSize);
    m_psiField.resize(m_gridSize + 1);

    for(int i = 0; i < m_gridSize; i++){
        m_sigmaField[i].resize(m_gridSize + 1);
        m_sigmaField[i].fill(INIT_SIGMA_GAMMA);

        m_gammaField[i].resize(m_gridSize);
        m_gammaField[i].fill(INIT_SIGMA_GAMMA);

        m_phiField[i].resize(m_gridSize + 1);
        m_phiField[i].fill(INIT_PHI_PSI);

        m_psiField[i].resize(m_gridSize);
        m_psiField[i].fill(INIT_PHI_PSI);
    }
<<<<<<< HEAD:final_qt/water/fluidGPU.cpp
*/
    //initialize sigma and gamma inside the dampening region
/*    for(int i = 0; i < m_gridSize; i++){
        for(int j= 0; j < m_gridSize; j++){
=======
    m_gammaField[m_gridSize].resize(m_gridSize);
    m_gammaField[m_gridSize].fill(INIT_SIGMA_GAMMA);
    m_psiField[m_gridSize].resize(m_gridSize);
    m_psiField[m_gridSize].fill(INIT_PHI_PSI);

    //initialize sigma and gamma inside the dampening region
    for(int i = 0; i <= m_gridSize; i++){
        for(int j = 0; j <= m_gridSize; j++){
>>>>>>> 45d010989ed255ec99ace55830acd9541d0b8402:final_qt/water/fluid.cpp
            if(i < DAMPENING_REGION || i >= m_gridSize - DAMPENING_REGION ||
                    j < DAMPENING_REGION || j >= m_gridSize - DAMPENING_REGION){
                //horizontal and vertical distances
                float iDistance = 0;
                float jDistance = 0;
                if(i < DAMPENING_REGION){
                    iDistance = DAMPENING_REGION - i;
                } else if(i >= m_gridSize - DAMPENING_REGION){
                    iDistance = i - (m_gridSize - DAMPENING_REGION - 1);
                }

                if(j < DAMPENING_REGION){
                    jDistance = DAMPENING_REGION - j;
                } else if(j >= m_gridSize - DAMPENING_REGION){
                    jDistance = j - (m_gridSize - DAMPENING_REGION - 1);
                }

                iDistance /= (float)DAMPENING_REGION;
                jDistance /= (float)DAMPENING_REGION;

                //distance
                float distance = sqrt((iDistance * iDistance) + (jDistance * jDistance));

                //quadratic function
                float value = (QUADRATIC_A * distance * distance) + (QUADRATIC_B * distance) + QUADRATIC_C;

<<<<<<< HEAD:final_qt/water/fluidGPU.cpp
                const int index = i*m_gridSize + j;
                m_sigmaField[index] = value;
                m_gammaField[index] = value;
=======
                if(i < m_gridSize){
                    m_sigmaField[i][j] = value;
                }
                if(j < m_gridSize){
                    m_gammaField[i][j] = value;
                }
>>>>>>> 45d010989ed255ec99ace55830acd9541d0b8402:final_qt/water/fluid.cpp
            }
        }
    }
*/
#ifdef Fluid_DEBUG
   // m_gridSize = 20;
    const float END_TIME = 1;
    m_dt = 0.2;
    m_dx = 0.5;
    m_dxInv = 1/0.5;
    int count = 0;
    float timeCount = 0.f;
    while( timeCount < END_TIME )
    {
        update(m_dt);
        if( m_updateCount %20 == 0 )
            writeToImage(VELOCITY);
        count++;
        printf("Count:%d\n", count);
       printMat( m_depthField );
        timeCount+= m_dt;
    }
    assert(0);
#endif

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
        glPushMatrix();
        glBegin(GL_POINTS);
        glColor4f( m_color.r,m_color.g, m_color.b, m_color.a );
        for( int i = 0; i < m_gridSize; i++ )
        {
            for( int j =0; j < m_gridSize; j++ )
            {
                float posX = -m_domainSize+ j;
                float posZ = -m_domainSize + i;
                glVertex3f(posX, m_heightField[getIndex1D(i,j,HEIGHT)],posZ);
            }
        }
        glPopMatrix();
        glEnd();
    }
    else if ( method == DRAW_MESH )
    {
        // This mode doesn't deal with the invisible triangles
        glPushMatrix();
        glEnable(GL_BLEND);
        glBlendFunc(GL_ONE,GL_ONE);
        float halfDomain = m_domainSize/2;

        const bool drawStrip = false;
        if( drawStrip == true )
        {
            for( int i = 0; i < m_gridSize-1; i++ )
            {

                glBegin( GL_TRIANGLE_STRIP );
                glColor4f(m_color.r,m_color.g,m_color.b,m_color.a);
                for( int j = 0; j < m_gridSize; j++ )
                {
                    int index = i*m_gridSize + j;
                    glNormal3f( m_normalField[index].x, m_normalField[index].y, m_normalField[index].z );
                    glVertex3f(  -halfDomain+j*m_dx, m_heightField[index], -halfDomain+i*m_dx );
                    index = (i+1)*m_gridSize + j;
                    glNormal3f( m_normalField[index].x, m_normalField[index].y, m_normalField[index].z );
                    glVertex3f(  -halfDomain +j*m_dx, m_heightField[index] , -halfDomain +(i+1)*m_dx );
                }
                glEnd();
            }
        }
        else
        {
             // Else we use triangles mode to draw(This will hide the invisible triangles)
            glBegin( GL_TRIANGLES );
            glColor4f(m_color.r,m_color.g,m_color.b,m_color.a);
            for( int i = 0;  i < m_triangles.size(); i++ )
            {
                Tri curTri  = m_triangles[i];
                int r[3]; int c[3];
                r[0] = curTri.a2D.indRow; c[0] = curTri.a2D.indCol;
                r[1] = curTri.b2D.indRow; c[1] = curTri.b2D.indCol;
                r[2] = curTri.c2D.indRow; c[2] = curTri.c2D.indCol;
                // Firstly check if the triangle is visible
                int count = 0;

                for( int m = 0; m < 3; m++ )
                {
                    if( m_depthField[getIndex1D(r[m],c[m],DEPTH)] > EPSILON )
                        count++;
                }
                if( count == 0 )
                    continue;

                float tx, ty, tz;
                for( int m = 0; m < 3; m++ )
                {
                    //if( m_depthField[getIndex1D(r[m],c[m],DEPTH)] > EPSILON )
                      glColor4f(m_color.r,m_color.g,m_color.b,m_color.a);
                   /*else
                        glColor4f(0.f,0.f,0.f, m_color.a );
                    */
                    int index = r[m]*m_gridSize + c[m];
                    tx = -halfDomain + c[m]*m_dx;
                    ty = m_heightField[index];
                    tz = - halfDomain + r[m]*m_dx;
                    glNormal3f( m_normalField[index].x, m_normalField[index].y, m_normalField[index].z );
                    glVertex3f( tx, ty, tz );
                }
            }
            glEnd();
        }
        glDisable(GL_BLEND);
        glPopMatrix();
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

        const float halfDomain = m_domainSize/2;
        const int magn = 2;
        for (int row = 0; row < m_gridSize; row++)
        {
            for (int column = 0; column < m_gridSize; column++)
            {
                glBegin(GL_LINES);

                Vector3 curVert = Vector3(-halfDomain+ column*m_dx,
                                          m_heightField[getIndex1D(row,column,HEIGHT)],
                                          -halfDomain +row*m_dx);
                Vector3 curNorm = m_normalField[getIndex1D(row,column,NORMAL)];

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
            Tri t1,t2;
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
    /*
    float hRest = computeHRest();

    //iterate through the dampening region
    for(int i = 1; i < m_gridSize - 1; i++){
        for(int j = 1; j < m_gridSize - 1; j++){
            if(i < DAMPENING_REGION || i >= m_gridSize - DAMPENING_REGION ||
                    j < DAMPENING_REGION || j >= m_gridSize - DAMPENING_REGION){
                // Equation 10
                // h(i,j) += ((-sigma(i,j) * (h(i,j) - hRest)) + phi(i,j)) * delta_t
                // Equation 21
                // h(i,j) += ((-gamma(i,j) * (h(i,j) - hRest)) + psi(i,j)) * delta_t
                float currH = m_depthField[i][j] + m_terrainHeightField[i][j];
                float eq10 = ((-m_sigmaField[i][j] * (currH - hRest)) + m_phiField[i][j]) * m_dt;
                float eq21 = ((-m_gammaField[i][j] * (currH - hRest)) + m_psiField[i][j]) * m_dt;
                m_depthField[i][j] += eq10 + eq21;

                // Equation 11
                // u(i+0.5,j) += -0.5 * (sigma(i+1,j) + sigma(i,j)) * u(i+0.5,j) * delta_t
                m_velocityU[i][j] += -0.5 * (m_sigmaField[i + 1][j] + m_sigmaField[i][j]) * m_velocityU[i][j] * m_dt;

                // Equation 22
                // w(i,j+0.5) += -0.5 * (gamma(i,j+1) + gamma(i,j)) * w(i,j+0.5) * delta_t
                m_velocityW[i][j] += -0.5 * (m_gammaField[i][j + 1] + m_gammaField[i][j]) * m_velocityW[i][j] * m_dt;

                // Equation 12
                // phi(i,j) += -LAMBDA_UPDATE * sigma(i,j) * ((w(i,j+0.5) - w(i,j-0.5)) / delta_x) * delta_t
                m_phiField[i][j] += -LAMBDA_UPDATE * m_sigmaField[i][j] * (m_velocityW[i][j] - m_velocityW[i][j - 1]) *
                        m_dxInv * m_dt;

                // Equation 13
                // phi(i,j) *= LAMBDA_DECAY
                m_phiField[i][j] *= LAMBDA_DECAY;

                // Equation 23
                // psi(i,j) += -LAMBDA_UPDATE * psi(i,j) * ((u(i+0.5,j) - u(i-0.5,j)) / delta_x) * delta_t
                // CORRECTION: psi(i,j) should be gamma(i,j) on the right-hand side
                m_psiField[i][j] += -LAMBDA_UPDATE * m_gammaField[i][j] * (m_velocityU[i][j] - m_velocityU[i - 1][j]) *
                        m_dxInv * m_dt;

                // Equation 24
                // psi(i,j) *= LAMBDA_DECAY
                m_psiField[i][j] *= LAMBDA_DECAY;
            }
        }
    }*/
}

/**
 * @brief computes the resting height of the fluid
 * the resting height is computed as the average across the depth field
 * the result of this function is used in dampenWaves()
 * @return h_rest
 */
float FluidGPU::computeHRest(){
    float hRest = 0;
/*
    for(int i = 0; i < m_gridSize; i++){
        for(int j = 0; j < m_gridSize; j++){
            hRest += m_depthField[i][j] + m_terrainHeightField[i][j];
        }
    }
*/
    return (hRest / (float)(m_gridSize * m_gridSize));
}

void FluidGPU::clampFields(){
    /*
    // clamp h(i,j) >= 0
    // clamp u(i,j) < alpha * (delta_x / delta_t)
    // clamp w(i,j) < alpha * (delta_x / delta_t)

    float velocityClamp = CLAMP_ALPHA * (m_dx / m_dt);

    for(int i = 0; i < m_gridSize; i++){
        for(int j = 0; j < m_gridSize; j++){
            m_depthField[i][j] = max(0.0f, m_depthField[i][j]);
            m_velocityU[i][j] = min(velocityClamp, m_velocityU[i][j]);
            m_velocityW[i][j] = min(velocityClamp, m_velocityW[i][j]);
        }

        m_velocityU[i][m_gridSize] = min(velocityClamp, m_velocityU[i][m_gridSize]);
        m_velocityW[m_gridSize][i] = min(velocityClamp, m_velocityW[m_gridSize][i]);
    }
    */
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
    case HEIGHT:
    case NORMAL:
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
    copybackGPU(DEPTH,m_depthField);
    copybackGPU(HEIGHT,m_heightField);
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
    glBegin(GL_POINTS);
    //glEnable(GL_PROGRAM_POINT_SIZE_EXT);
    //glPointSize(20);
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

    for(int i = 0; i < NUM_DROPPING_PARTICLES; i++){
        //jitter particle positions
        float randX = randomFloatGenerator(-radius, radius);
        float randY = randomFloatGenerator(0, PARTICLE_DROP_RANGE);
        float randZ = randomFloatGenerator(-radius, radius);

        Vector3 position = Vector3(coordX + randX, PARTICLE_DROP_HEIGHT + randY, coordY + randZ);
        Vector3 velocity = Vector3::zero();
        Vector3 acceleration = Vector3(0, GRAVITY, 0);

        //make new particle
        Particle *newParticle = new Particle(SPLASH_PARTICLE_RADIUS, Veff, position, velocity, acceleration);
        m_particles.append(newParticle);
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
