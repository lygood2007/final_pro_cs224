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
    safeFreeArray1D( m_tempBuffer );
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

    advect( DEPTH, m_depthField );
    advect( VEL_U, m_velocityU );
    advect( VEL_W, m_velocityW );

    updateDepth();
    updateHeight();
    updateVelocities();

    applyBoundary();
    checkBoundary();

    //clampFields();

#ifdef DAMPEN_WAVES
    //dampenWaves();
#endif

    computeNormal();

    m_timeElapsed += dt;
    m_updateCount++;

#ifdef SAVE_IMAGE
    const int savePerFrames = 5;
    if( (m_updateCount%savePerFrames) == 0 )
        saveToImage( VELOCITY );
#endif

}

/**
 * @brief Increment the height of the rectangular region around (posX, posZ) by incHeight
 * @param posX The x position
 * @param posZ The z position
 * @param radius The radius
 * @param incHeight The amount of height added
 */
void FluidGPU::incrementH(const int posX, const int posZ, const int radius, const float incHeight )
{
    for( int i = posZ - radius; i < posZ + radius + 1; i++ )
    {
        for( int j = posX - radius; j < posX + radius + 1; j++ )
        {
            if( i < 0 || i >= m_gridSize  || j < 0 || j >= m_gridSize )
            {
                continue;
            }
            else
            {
                const int index = getIndex1D( i,j,DEPTH );
                m_depthField[index] += incHeight;
                if( m_depthField[index] > maxHeight )
                    m_depthField[index] = maxHeight;
            }
        }
    }
}

void FluidGPU::addDrop(const int posX, const int posZ)
{
    // Fixed size
    int radius = m_gridSize/20;
    if( radius < 1 )
        radius = 1;
    float incH = 8;
    incrementH( posX, posZ, radius, incH );
}

/**
 * @brief Add drop to random positions
 *
 */
void FluidGPU::addRandomDrop( const float freq )
{
    // freq is not useful right now
    //float rnd = randomFloatGenerator();
        //int posX = rand()%m_gridSize;
        //int posY = rand()%m_gridSize;
       // int radius = rand()%(m_gridSize/15);
        /**
         * @brief not applicable right now
         */
        int posX = 0.5*m_gridSize;
        int posZ = 0.5*m_gridSize;
        int radius = rand()%(m_gridSize/15);
    float rndHeight = randomFloatGenerator( 5, 15 );

    incrementH( posX, posZ, radius, rndHeight );
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
    assert( t->getBound() == m_domainSize/2 );
   const Vector3* terrainVertices = t->getVerts();
   for( int i = 0; i < m_gridSize; i++ )
   {
       for( int j = 0; j < m_gridSize; j++ )
       {
           const int index = getIndex1D( i,j, T_HEIGHT);
           m_terrainHeightField[index] = terrainVertices[index].y;
       }
   }
    initDepthField();
    updateHeight();
//   printMat( m_terrainHeightField );
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

    size = (m_gridSize+1)*(m_gridSize+1)*sizeof(float);
    m_tempBuffer = (float*)malloc(size);
    for( int i = 0; i < (m_gridSize+1)*(m_gridSize+1); i++ )
    {
        m_tempBuffer[i] = 0.f;
    }
    size = m_gridSize*m_gridSize*sizeof(Vector3);
    m_normalField = (Vector3*)malloc(size);
    for( int i = 0; i < m_gridSize*m_gridSize; i++ )
    {
        // Initial value for normals
        m_normalField[i] = Vector3(0,1,0);
    }

    /*
    m_depthField.resize(m_gridSize);
    m_normalField.resize( m_gridSize );
    m_terrainHeightField.resize( m_gridSize );
    m_velocityU.resize(m_gridSize);
    m_velocityW.resize(m_gridSize+1);

    for( int i = 0; i < m_gridSize; i++ )
    {
        m_depthField[i].resize(m_gridSize);
        //m_depthField[i].fill(defaultHeight);
        m_normalField[i].resize(m_gridSize);
        m_normalField[i].fill( Vector3(0.f,1.f,0.f));

        m_terrainHeightField[i].resize(m_gridSize);
        m_terrainHeightField[i].fill(0.0);
        m_velocityU[i].resize(m_gridSize+1);
        m_velocityU[i].fill(defaultU);

        m_velocityW[i].resize(m_gridSize);
        m_velocityW[i].fill(defaultW);
    }
    m_velocityW[m_gridSize].resize(m_gridSize);
    m_velocityW[m_gridSize].fill(defaultW);*/

/*    assert( m_depthField.size() == m_gridSize && m_depthField[0].size() == m_gridSize );
    assert( m_velocityU.size() == m_gridSize && m_velocityU[0].size() == m_gridSize+1 );
    assert( m_velocityW.size() == m_gridSize+1&& m_velocityW[0].size() == m_gridSize );
*/

    buildTriangleList();

/*
    m_tempBuffer = new float*[m_gridSize+1];
    for( int i = 0; i < m_gridSize + 1; i++ )
    {
        m_tempBuffer[i] = new float[m_gridSize+1];
    }
*/
    // Set the random seed
    srand((unsigned)time(0));

    m_updateCount = 0;
    m_timeElapsed = 0.f;
    // Default color
    m_color = Colorf(0.1f,0.4f,0.8f,1.0f);
    m_renderNormals = false;

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

#ifdef USE_CUDA

#endif
}

/**
 * @brief Advect the fluid
 */
void FluidGPU::advect( FieldType type, float* vec  )
{
    int width;
    int height;
    switch( type )
    {
    case DEPTH:
        {
            width = m_gridSize;
            height = m_gridSize;
            break;
        }
     case VEL_U:
    {
        width = m_uWidth;
        height = m_gridSize;
        break;
    }
    case VEL_W:
    {
        width = m_gridSize;
        height = m_gridSize+1;
        break;
    }
    default:
        assert(0);
        break;
    }

    for( int i = 1; i < height - 1; i++ )
    {
        for( int j = 1; j < width - 1; j++ )
        {
            float u,w;
            switch( type )
            {
            case DEPTH:
            {
                u = 0.5*( m_velocityU[getIndex1D(i,j,VEL_U)]+m_velocityU[getIndex1D(i,j+1,VEL_U)] );
                w = 0.5*( m_velocityW[getIndex1D(i,j,VEL_W)]+m_velocityW[getIndex1D(i+1,j,VEL_W)] );
                break;
            }
            case VEL_U:
            {
                u = m_velocityU[getIndex1D(i,j,VEL_U)];
                w = 0.25*( m_velocityW[getIndex1D(i,j,VEL_W)]+m_velocityW[getIndex1D(i,j-1,VEL_W)]
                           +m_velocityW[getIndex1D(i+1,j-1,VEL_W)]+m_velocityW[getIndex1D(i+1,j,VEL_W)] );
                break;
            }
            case VEL_W:
            {
                u = 0.25*( m_velocityU[getIndex1D(i-1,j,VEL_U)]+m_velocityU[getIndex1D(i-1,j+1,VEL_U)]
                           +m_velocityU[getIndex1D(i,j,VEL_U)]+m_velocityU[getIndex1D(i,j+1,VEL_U)] );
                w = m_velocityW[getIndex1D(i,j,VEL_W)];
                break;
            }
            default:
            {
                assert(0);
                break;
            }
            }
            float curPosX = (float)j;
            float curPosY = (float)i;
            float prev_x = curPosX - u*m_dt*m_dxInv;
            float prev_z = curPosY - w*m_dt*m_dxInv;
           m_tempBuffer[getIndex1D(i,j,TMP)] = bilinearInterp( vec, prev_x, prev_z, width, height );
        }
    }

    // Copy back
    for( int i = 1; i < height - 1; i++ )
    {
        for( int j = 1; j < width - 1; j++ )
        {
            switch( type )
            {
            case DEPTH:
            {
                vec[getIndex1D(i,j,DEPTH)] = m_tempBuffer[getIndex1D(i,j,TMP)];
                break;
            }
            case VEL_W:
            {
                vec[getIndex1D(i,j,VEL_W)] = m_tempBuffer[getIndex1D(i,j,TMP)];
                break;
            }
            case VEL_U:
            {
                vec[getIndex1D(i,j,VEL_U)] = m_tempBuffer[getIndex1D(i,j,TMP)];
                break;
            }
            default:
                assert(0);
            }
        }
    }
}

/**
 * @brief updateVelocities Update the velocities
 */
void FluidGPU::updateVelocities()
{
    float h1,h2;
    for (int i=1;i<m_gridSize-1;i++)
    {
        for (int j=2;j<m_gridSize-1;j++)
        {
            //h1 = max(0.0f, m_terrainHeightField[getIndex1D(i,j,T_HEIGHT)] + m_depthField[getIndex1D(i,j,DEPTH)] );
            //h2 = max(0.0f, m_terrainHeightField[getIndex1D(i,j-1,T_HEIGHT)] + m_depthField[getIndex1D(i,j-1,T_HEIGHT)] );
            h1 = m_heightField[getIndex1D(i,j,HEIGHT)];
            h2 = m_heightField[getIndex1D(i,j-1,HEIGHT)];
            m_velocityU[getIndex1D(i,j,VEL_U)] += GRAVITY * m_dt * m_dxInv * ( (h1-h2) );
        }
    }
    for (int i=2;i<m_gridSize-1;i++)
    {
        for (int j=1;j<m_gridSize-1;j++)
        {
            h1 = m_heightField[getIndex1D(i,j,HEIGHT)];
            h2 = m_heightField[getIndex1D(i-1,j,HEIGHT)];
            //h1 = max(0.0f, m_terrainHeightField[getIndex1D(i,j,T_HEIGHT)] + m_depthField[getIndex1D(i,j,DEPTH)] );
            //h2 = max(0.0f, m_terrainHeightField[getIndex1D(i-1,j,T_HEIGHT)] + m_depthField[getIndex1D(i-1,j,DEPTH)] );
            m_velocityW[getIndex1D(i,j,VEL_W)] += GRAVITY* m_dt * m_dxInv * ((h1-h2));
        }
    }
}

/**
 * @brief updateDepth Update the depth field
 */
void FluidGPU::updateDepth()
{
    const float decay = 1.f;

    for( int i = 1; i < m_gridSize - 1; i++ )
    {
        for( int j = 1; j < m_gridSize - 1; j++ )
        {
            float dh = -decay*m_depthField[getIndex1D(i,j,DEPTH)]*m_dxInv*((m_velocityU[getIndex1D(i,j+1,VEL_U)] - m_velocityU[getIndex1D(i,j,VEL_U)])
                                                           +(m_velocityW[getIndex1D(i+1,j,VEL_W)] - m_velocityW[getIndex1D(i,j,VEL_W)]));
            m_depthField[getIndex1D(i,j,DEPTH)] += dh*m_dt;
        }
    }
}

/**
 * @brief updateHeight Update the height field
 */
void FluidGPU::updateHeight()
{
    for( int i = 0; i < m_gridSize*m_gridSize; i++ )
    {
        m_heightField[i] = max(0.f, m_depthField[i] + m_terrainHeightField[i]);
    }
}

/**
 * @brief Apply boundary condition
 */
void FluidGPU::applyBoundary()
{
    for( int i = 0; i < m_gridSize; i++ )
    {
        m_depthField[getIndex1D(0,i,DEPTH)] = max(0.f, m_heightField[getIndex1D(1,i,HEIGHT)]
                                                   - m_terrainHeightField[getIndex1D(0,i,T_HEIGHT)]);
        m_depthField[getIndex1D(m_gridSize-1,i,DEPTH)] = max( 0.f, m_heightField[getIndex1D(m_gridSize - 2,i,HEIGHT)]
                                                              - m_terrainHeightField[getIndex1D(m_gridSize-1,i,T_HEIGHT)]);
    }

    for( int j = 0; j < m_gridSize; j++ )
    {
        m_depthField[getIndex1D(j,0,DEPTH)] = max(0.f, m_heightField[getIndex1D(j,1,HEIGHT)] - m_terrainHeightField[getIndex1D(j,0,T_HEIGHT)]);
        m_depthField[getIndex1D(j,m_gridSize-1,DEPTH)] = max(0.f, m_heightField[getIndex1D(j,m_gridSize-2,HEIGHT)]
                                                            - m_terrainHeightField[getIndex1D(j,m_gridSize-1,T_HEIGHT)]);
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
 * @brief Computet the normal of each points
 */
void FluidGPU::computeNormal()
{
    for( int i = 0; i < m_gridSize; i++ )
    {
        for( int j = 0; j < m_gridSize; j++ )
        {   
            int numNeighbors = 0;
            Vector3 offsets[8];
            // Search fir eight neightbors
            Vector2 coords[8];
            QList<Vector2 > neighbors;
            coords[0] = Vector2(i,     j - 1);
            coords[1] = Vector2(i + 1, j - 1);
            coords[2] = Vector2(i + 1, j);
            coords[3] = Vector2(i + 1, j + 1);
            coords[4] = Vector2(i,     j + 1);
            coords[5] = Vector2(i - 1, j + 1);
            coords[6] = Vector2(i - 1, j);
            coords[7] = Vector2(i - 1, j - 1);
            for( int m = 0; m < 8; m++ )
            {
                if( coords[m].x < 0 || coords[m].y < 0 || coords[m].x > m_gridSize - 1 || coords[m].y > m_gridSize - 1 )
                    continue;
                neighbors.push_back( coords[m]);
            }

            for( int m = 0; m < neighbors.size(); m++ )
            {
                //offsets[m] = Vector3(neighbors[m].x,m_depthField[neighbors[m].x][neighbors[m].y],neighbors[m].y) - Vector3(i,m_depthField[i][j],i);
                offsets[m] = Vector3(neighbors[m].x,m_heightField[getIndex1D(neighbors[m].x,neighbors[m].y,HEIGHT)],neighbors[m].y) -
                        Vector3(i,m_heightField[getIndex1D(i,j,HEIGHT)]
                                ,j);
            }

            Vector3 sum = Vector3::zero();
            for (int m = 0; m < neighbors.size(); ++m)
            {
                Vector3 tmp = Vector3::zero();
                if( m+1 == numNeighbors )
                    tmp = offsets[m].cross(offsets[0]);
                else
                    tmp = offsets[m].cross(offsets[m+1]);

                sum += tmp;
            }

            m_normalField[getIndex1D(i,j,NORMAL)] = -sum.getNormalized();
        }
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
        const int magn = 5;
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
 * @brief initDepthField Initialize the depth field
 */
void FluidGPU::initDepthField()
{
    for( int i = 0; i < m_gridSize; i++ )
    {
        for( int j =0; j < m_gridSize; j++ )
        {
            m_depthField[getIndex1D(i,j,DEPTH)] = max(0.f, defaultHeight - m_terrainHeightField[getIndex1D(i,j,T_HEIGHT)] );
        }
    }
}

/**
 * @brief build the triangle List
 */
void FluidGPU::buildTriangleList()
{
    /*
    assert( m_depthField.size() == m_gridSize );
    assert( m_depthField[0].size() == m_gridSize );
    */

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
    case T_HEIGHT:
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
