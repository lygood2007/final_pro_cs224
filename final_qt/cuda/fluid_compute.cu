/** fluid_compute.cu
 ** Brief: Deal with all the computation here
 ** Project: large-scale fluids
 ** Date: 04/10/2013
 ** Member: Scott, Hobarts, Yan Li
 **/

#ifndef FLUID_COMPUTE_CU
#define FLUID_COMPUTE_CU

#include <cuda.h>
#include <stdio.h>
#include <assert.h>
#include "fluid_global.h"
// includes, cuda
#include <cuda_runtime.h>
#include <cuda_gl_interop.h>

#define EPS 0.001
extern "C"
{
void initGridGPU( const int hostGridSize, const float hostdx, const float halfdm, const float* hostTerrainMap );
void copybackGPU(FieldType type, float* hostMap  );
void destroyGPUmem();
void addDropGPU(const int posX, const int posZ, const int radius, const float h );
void advectGPU(const float dt);
void updateFluidGPU( const float dt );
bool findSupportDevice();

void initParticlesGPU(const int numParticles);
void updateParticlesGPU( const float dt, const float accX, const float accY, const float accZ );
void intersectParticlesGPU( const float halfDomain, const float mdxInv, const float heightChange );
void inputParticlesGPU(const float *particlePositions, const float *particleVelocities );

void clampFieldsGPU( const float velocityClamp );

void initDampeningFieldsGPU( const int sizeDampeningRegion, const float quadraticA, const float quadraticB, const float quadraticC );
void dampenWavesGPU( const float hRest, const float dt, const float dxInv, const float lambdaUpdate, const float lambdaDecay );
}

__host__ __device__ inline float cudaMax( float a, float b )
{
    return (a>b)?(a):(b);
}

__host__ __device__ inline float cudaMin( float a, float b )
{
    return (a<b)?(a):(b);
}


// The 2D vector structure for GPU computing
 struct vec2
{
    union
    {
        struct {float x,y; };
        float xy[2];
    };
};

 struct vec3
{
    union
    {
        struct {float x,y,z; };
        float xyz[3];
    };
};

const int blockSizeX = 4;
const int blockSizeY = 4;

vec3* devicePaintMap; // Paint map, we copy back this buffer for drawing. It stores the position of vertices
vec3* deviceNormalMap; // Normal map for GPU
float* deviceTerrainMap; // Terrain map for GPU
float* deviceHeightMap; // Height map for GPU
float* deviceDepthMap; // Depth map for GPU
float* deviceVelocityUMap; // VelocityU map for GPU
float* deviceVelocityWMap; // VelocityW map for GPU

float* deviceNextDepthMap; // Temp buffer for storing next depth map
float* deviceNextVelocityUMap; // Temp buffer for storing next velocity U map
float* deviceNextVelocityWMap; // Temp buffer for storing next velocity W map

//particle data structures
vec3* deviceParticlePositionsArray; // particle positions array
vec3* deviceParticleVelocitiesArray; // particle velocities array
int deviceNumParticles; // number of particles

//dampening waves data structures
float* deviceSigmaMap;
float* deviceGammaMap;
float* devicePhiMap;
float* devicePsiMap;
int deviceDampeningRegion; // size of dampening region

/**
 * pitches for the maps above
 */
// Error
cudaError_t error;

// The grid size for heightmap, depthmap, terrainmap
int gridSize;
// The width for velocity u
int uwidth;
// The height for veloctiy u
int uheight;
// The width for velocity w
int wwidth;
// The height for veloctiy w
int wheight;

// dx
float mapdx;
// dxInv
float mapdxInv;
// halfDomain
float halfDomain;

// Forward decaration
void initGridGPU( const int hostGirdSize, const float* hostTerrainMap );
void addDropGPU(const int posX, const int posZ, const int radius, const float h );
void copybackGPU(float* hostHeightMap );
void destroyGPU();
void advectGPU(const float dt);
void updateFluidGPU( const float dt );

//particles forward declarations
void initParticlesGPU(const int numParticles);
void updateParticlesGPU( const float dt, const float accX, const float accY, const float accZ );
void intersectParticlesGPU( const float halfDomain, const float mdxInv, const float heightChange );
void inputParticlesGPU( const float *particlePositions, const float *particleVelocities );

//stabilize forward declarations
void clampFieldsGPU( const float velocityClamp );

//dampening waves forward declarations
void initDampeningFieldsGPU( const int sizeDampeningRegion, const float quadraticA, const float quadraticB, const float quadraticC );
void dampenWavesGPU( const float hRest, const float dt, const float dxInv, const float lambdaUpdate, const float lambdaDecay );

void checkInitializedDeviceField( float* device, int width, int height )
{
    float* host = (float*)malloc(width*height*sizeof(float));
    cudaMemcpy( host, device, width*height*sizeof(float),cudaMemcpyDeviceToHost);

    for( int i = 0; i < width*height; i++ )
    {
        if( host[i] != 0.f )
        {
            assert(0);
        }
    }
    free(host);
}

void  checkCudaError( cudaError_t error )
{
    if( error != cudaSuccess )
    {
        //cout <<"CUDA error code: "<<cudaGetErrorString(error);
        printf( "CUDA error code: %s\n",cudaGetErrorString(error) );
    }
}
/**
 * Check if the pointer is null, if it's null, exit the program
 */
template <class T>
void check1DNotNull(T* array )
{
    if( array == NULL )
    {
        printf("Wrong pointer\n");
    }
}

/**
 * Initialize a vec2
 */
__host__ __device__ inline vec2 initVec2( float x, float y )
{
    vec2 result;
    result.x = x;result.y = y;
    return result;
}

/**
 *  Initialize a vec3
 */
__host__ __device__ inline vec3 initVec3( float x, float y, float z )
{
    vec3 result;
    result.x = x;result.y = y; result.z = z;
    return result;
}

// Review passed
/**
* Compute the cross product of two vectors
*/
__host__ __device__ inline vec3 cross( const vec3 v1, const vec3 v2 )
{
    vec3 result = initVec3(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
    return result;
}

// Review passed
/**
 * Compute the normalized vector
 */
__host__ __device__ inline vec3 normalize( const vec3 v )
{
    float d = sqrt(v.x*v.x + v.y*v.y + v.z*v.z );
    return initVec3( v.x/d,v.y/d, v.z/d );
}

// Review passed
/**
 * Functions for getting data from 2D array in GPU
 */
__host__ __device__ inline float map2Dread( const float* map, int i, int j, int width )
{
    return map[i*width+ j];
}
// Review passed
/**
 * Functions for writting data from 2D array in GPU
 */
__host__ __device__ inline void map2Dwrite( float* map, int i, int j, float value, int width )
{
    map[i*width + j] = value;
}

// Review passed
/**
 *  Initialize the depth
 */
__global__ void initDepthCUDA( float* depthMap, const float* terrainMap, const int width, const int height )
{
    int i = blockDim.y*blockIdx.y +threadIdx.y;
    int j = blockDim.x*blockIdx.x +threadIdx.x;

    if( i >= 0 && i < height && j >= 0 && j < width )
    {
        float depth = cudaMax(0.f,defaultHeight - map2Dread( terrainMap, i,j, width ));
       // float depth = 5;
        map2Dwrite( depthMap,i,j,depth,width );
    }
}

// Review passed
/**
 * Update the height field by plusing the depth and terrain
 */
__global__ void updateHeightCUDA( float*heightMap, const float* depthMap, const float* terrainMap,
                              const int width, const int height )
{
    int i = blockDim.y*blockIdx.y + threadIdx.y;
    int j = blockDim.x*blockIdx.x + threadIdx.x;

    if( i >= 0 && i < height && j >= 0 && j < width )
    {
        float depth = map2Dread( depthMap,i,j,width );
        float terrainHeight = map2Dread( terrainMap, i,j,width );
        float h = depth + terrainHeight;
        map2Dwrite( heightMap, i,j, h ,width);
    }
}

// Review passed
/**
 * Add drop to specified rectangular region
 */
__global__ void addDropCUDA( float* depthMap, const int posX, const int posZ, const int radius,
                             const float h, const int width, const int height )
{
    int i = blockDim.y*blockIdx.y + threadIdx.y;
    int j = blockDim.x*blockIdx.x + threadIdx.x;

    if( i>= cudaMax(posZ-radius,0) && i < cudaMin(posZ+radius+1,height)
            && j >= cudaMax(posX-radius,0)&&j < cudaMin(posX+radius+1,width)
            )
    {
        float newH = cudaMin(map2Dread(depthMap,i,j,width)+h,maxHeight );
        map2Dwrite( depthMap, i,j, newH, width );
    }
}

// Review passed
/**
 * bilinear interpolation
 */
 __host__ __device__ float bilinearIerp( const float* vec, float x, float z, const int width, const int height )
 {
     if( x < 0 )
         x = 0.f;
     if( z < 0 )
         z = 0.f;
     if( x > width - 1 )
         x = width - 1;
     if( z > height - 1 )
         z = height -1;

     const int X = (int)x;
     const int Y = (int)z;
     const float s1 = x - X;
     const float s0 = 1.f - s1;
     const float t1 = z - Y;
     const float t0 = 1.f-t1;
     float e1, e2, e3,e4;
     e1 = e2 = e3 = e4 = 0;
     //e1 = vec[Y*width+X];
     e1 = map2Dread(vec,Y,X,width);
     if( Y+1 <= height- 1 )
     {
      //   e2 = vec[(Y+1)*width + X];
         e2 = map2Dread( vec, Y+1,X,width);
     }
     if( X +1 <= width -1 )
     {
        // e3 = vec[Y*width + X+1];
         e3 = map2Dread( vec, Y, X+1, width );
     }
     if( Y+1 <= height - 1 && X + 1 <= width - 1)
     {
      //   e4 = vec[(Y+1)*width + X+1];
         e4 = map2Dread( vec, Y+1,x+1,width );
     }

     float result = s0*(t0*e1 + t1*e2 )+
             s1*(t0*e3  + t1*e4 );

     return  result;
 }

 // Review passed
 /**
  * Advection: depth
  **/
 __global__ void advectDepthCUDA( const float* depthMap, float* nextDepthMap, const float* velUMap, const float* velWMap,
                              const int width, const int height, const float dt, const float dxInv)
{
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 1 && i < height -1 && j >= 1 && j < width - 1 )
     {
         // Be careful about the width of velU
         float uw = width + 1;
         float u = 0.5*(map2Dread(velUMap,i,j,uw) + map2Dread( velUMap,i,j+1,uw) );
         float w = 0.5*(map2Dread(velWMap,i,j,width) + map2Dread( velWMap,i+1,j,width) );

         float curPosX = (float)j;
         float curPosY = (float)i;
         float prev_x = curPosX - u*dt*dxInv;
         float prev_z = curPosY - w*dt*dxInv;
         map2Dwrite(nextDepthMap, i,j, bilinearIerp( depthMap,prev_x,prev_z,width,height ), width);
     }
 }

 // Review passed
 /**
  * Advection: velocity U
  */
 __global__ void advectVelUCUDA( const float* velUMap, float* nextVelUMap, const float* velWMap,
                            const int width, const int height, const float dt, const float dxInv )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 1 && i < height - 1&&j >= 1 && j < width - 1)
     {
         float ww = width - 1;
        float u = map2Dread( velUMap, i,j,width );
        float w = 0.25*(map2Dread(velWMap,i,j, ww ) + map2Dread( velWMap, i,j-1,ww) + map2Dread( velWMap,i+1,j-1,ww) + map2Dread( velWMap, i+1,j,ww ) );
        float curPosX = (float)j;
        float curPosY = (float)i;
        float prev_x = curPosX - u*dt*dxInv;
        float prev_z = curPosY - w*dt*dxInv;
        map2Dwrite(nextVelUMap, i,j, bilinearIerp(  velUMap,prev_x,prev_z,width,height ), width);
     }
 }

 // Review passed
 /**
  *  Advection: velocity W
  */
 __global__ void advectVelWCUDA( const float* velWMap, float* nextVelWMap, const float* velUMap,
                            const int width, const int height, const float dt, const float dxInv )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 1 && i < height - 1&&j >= 1 && j < width - 1)
     {
        float uw = width + 1;
        float u = 0.25*(map2Dread( velUMap,i,j,uw) + map2Dread(velUMap,i,j+1,uw) + map2Dread(velUMap,i-1,j+1,uw) + map2Dread(velUMap, i-1,j,uw) );
        float w = map2Dread( velWMap,i,j,width );
        float curPosX = (float)j;
        float curPosY = (float)i;
        float prev_x = curPosX - u*dt*dxInv;
        float prev_z = curPosY - w*dt*dxInv;
        map2Dwrite(nextVelWMap,i,j,bilinearIerp(  velWMap,prev_x,prev_z,width,height ), width);
     }
 }

 // Review passed
 /**
  * Update the depth field
  **/
 __global__ void updateDepthCUDA(float* depthMap, const float* velUMap, const float* velWMap,
                                 const int width, const int height, const float dt, const float dxInv )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 1 && i < height-1 && j >= 1 && j < width-1)
     {
            float decay = 1.f;
            float dep = map2Dread(depthMap,i,j,width);
            float dh = -decay*dep*dxInv*( (map2Dread(velUMap,i,j+1,width+1) - map2Dread(velUMap,i,j,width+1))
                                          + (map2Dread(velWMap,i+1,j,width) - map2Dread(velWMap,i,j,width)) );
            float nextDepth = dh*dt+dep;
            if( nextDepth < EPS )
            {
                 map2Dwrite(depthMap, i, j, 0.f, width );
            }
            else
            {
                map2Dwrite(depthMap, i, j, dh*dt+dep, width );
            }
     }
 }

 // Review passed
 /**
  * Update the velocity U field
  */
 __global__ void updateVelUCUDA( float* velUMap, const float* heightMap, const float* depthMap,
                                 const int width, const int height, const float dt, const float dxInv )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >=1 && i < height - 1 && j >= 2 && j < width - 2 )
     {
         // The width of heightmap is 1 smaller than the width of velocity U
        float h1 = map2Dread( heightMap, i,j, width - 1);
        float h2 = map2Dread(  heightMap, i,j-1, width-1 );
        float d1 = map2Dread(depthMap,i,j,width-1);
        float d2 = map2Dread(depthMap,i,j-1,width-1);

        // Read the origin value from velUMap
        float vel = map2Dread( velUMap, i,j, width );

        if( d1 < EPS || d2 < EPS )
        {
            float vel1 = map2Dread( velUMap,i,j-1,width);
            float vel2 = map2Dread( velUMap,i,j+1,width);
            float vel3 = map2Dread( velUMap, i,j, width );
            map2Dwrite( velUMap,i,j,0.33*(vel1+vel2+vel3),width);
            return;

        }
        float dv = GRAVITY*dt*dxInv*(h1-h2);

        // Add
        map2Dwrite( velUMap, i,j,vel+dv, width);
     }
     else
     {
         // for the bounday, we set the velocity to zero, using the neuman boundary condition
        map2Dwrite( velUMap, i,j,0.f, width);
     }
 }

 // Review passed
 /**
  * Update the velocity W field
  */
 __global__ void updateVelWCUDA( float* velWMap, const float* heightMap, const float* depthMap,
                                 const int width, const int height, const float dt, const float dxInv )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 2 && i < height - 2 && j >= 1 && j < width - 1 )
     {
         float h1 = map2Dread( heightMap, i,j, width );
         float h2 = map2Dread(  heightMap, i-1,j, width );

         float d1 = map2Dread(depthMap,i,j,width);
         float d2 = map2Dread(depthMap,i-1,j,width);

         float vel = map2Dread( velWMap, i,j, width );

         if( d1 < 0.0001 || d2 < 0.0001 )
         {
             float vel1 = map2Dread( velWMap,i-1,j,width);          
             float vel2 = map2Dread( velWMap,i+1,j,width);
             float vel3 = map2Dread( velWMap, i, j, width );
             map2Dwrite( velWMap,i,j,0.33*(vel1+vel2+vel3),width);
             return;

         }
         float dv = GRAVITY*dt*dxInv*(h1-h2);
         map2Dwrite( velWMap, i,j,vel+dv, width );
     }
     else
     {
         // for the bounday, we set the velocity to zero, using the neuman boundary condition
        map2Dwrite( velWMap, i,j,0.f, width);
     }
 }

 // Review passed
 /**
  * Apply the boundary condition
  **/
 __global__ void applyBoundaryCUDA( float* depthMap, const float* heightMap, const float* terrainMap,
                                    const int width, const int height )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
    float value;
    // Fix the boundary
    /*if( i == 0 || i == height-1|| j == 0 ||  j == width-1 )
    {
        map2Dwrite( depthMap, i,j, 0, width );
    }*/
     if( j == 0 && i !=  0 && i != height-1)
     {
         value = max(0.f, map2Dread( heightMap, i,1,width) - map2Dread( terrainMap, i,j,width ) );
         map2Dwrite( depthMap, i,j, value, width );
         return;
     }
     else if( j== width - 1&& i !=0 && i != height - 1 )
     {
         value = max( 0.f, map2Dread( heightMap, i, width - 2, width ) - map2Dread( terrainMap, i, j,width) );
        map2Dwrite( depthMap, i, j, value, width );
        return;
     }

     if( i == 0&& j != 0 && j != width -1 )
     {
         value = max(0.f, map2Dread( heightMap, 1, j, width) - map2Dread( terrainMap, i, j,width ));
         map2Dwrite( depthMap, i, j, value, width );
         return;
     }
     else if( i == height - 1&& j != width - 1 && j != 0)
     {
         value = max(0.f, map2Dread( heightMap, height - 2, j, width ) - map2Dread(terrainMap, i,j,width ));
         map2Dwrite( depthMap, i, j, value, width );
         return;
     }

     // Deal with the four courner, is there a way to simplify this? This function is too long!
     if( i== 0 && j == 0 )
     {
         value = max(0.f, map2Dread( heightMap, 1, 1, width ) - map2Dread(terrainMap, i,j,width ));
         map2Dwrite( depthMap, i,j, value, width );
         return;
     }
     else if( i==0 && j == width - 1 )
     {
         value = max(0.f, map2Dread( heightMap, 1, width-2, width ) - map2Dread(terrainMap, i,j,width ));
         map2Dwrite( depthMap, i,j, value, width );
         return;
     }
     else if ( i == height - 1 && j == 0 )
     {
         value = max(0.f, map2Dread( heightMap, height - 2, 1, width ) - map2Dread(terrainMap, i,j,width ));
         map2Dwrite( depthMap, i,j, value, width );
         return;
     }
     else if( i == height - 1 && j == width - 1 )
     {
         value = max(0.f, map2Dread( heightMap, height - 2, width-2, width ) - map2Dread(terrainMap, i,j,width ));
         map2Dwrite( depthMap, i,j, value, width );
         return;
     }
 }

 // Review passed
 /**
  * A general function to initialize the field to be all zero( Cannot use memset because we are using float )
  **/
 __global__ void initFieldCUDA( float* deviceMap, int width, int height )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 0&& i < height&&j >= 0 && j < width )
     {
         map2Dwrite( deviceMap, i, j, 0.f, width );
     }
 }

/**
 * Initialize the normal field
 **/
 __global__ void initNormalCUDA( vec3* normalMap, const int width, const int height )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 0&& i < height&&j >= 0 && j < width )
     {
         const int index = i*width + j;
         normalMap[index].x = 0;
         normalMap[index].y = 1;
         normalMap[index].z = 0;
     }
 }

 /**
  * Initialize the paint field
  */
 __global__ void updatePaintCUDA( vec3* paintMap, const float* heightMap, const float halfdm, const float dx, const int width, const int height )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 0&& i < height&&j >= 0 && j < width )
     {
         const int index = i*width + j;
         paintMap[index].x = -halfdm + j*dx;
         paintMap[index].y = heightMap[index];
         paintMap[index].z = -halfdm + i*dx;
     }
 }

 /**
  * compute the normals
  **/
 __global__ void computeNormalCUDA( vec3* normalMap, const float* heightMap, const int width, const int height )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 0&& i < height&&j >= 0 && j < width )
     {
            int numNeighbours = 0;

            vec3 offset[8];
            vec2 coords[8];
            vec2 neighbours[8];

            coords[0] = initVec2(i,     j - 1);
            coords[1] = initVec2(i + 1, j - 1);
            coords[2] = initVec2(i + 1, j);
            coords[3] = initVec2(i + 1, j + 1);
            coords[4] = initVec2(i,     j + 1);
            coords[5] = initVec2(i - 1, j + 1);
            coords[6] = initVec2(i - 1, j);
            coords[7] = initVec2(i - 1, j - 1);
            int m;
            for( m = 0; m < 8; m++ )
            {
                if( coords[m].x < 0 || coords[m].y < 0 || coords[m].x > height- 1 || coords[m].y > width- 1 )
                    continue;
                neighbours[numNeighbours] = coords[m];
                numNeighbours++;
            }

            for( m = 0; m < numNeighbours; m++ )
            {
                offset[m].x = neighbours[m].y - j;
                offset[m].z = neighbours[m].x - i;
                offset[m].y = map2Dread(heightMap,neighbours[m].x,neighbours[m].y,width) - map2Dread( heightMap, i,j,width );
            }

            vec3 sum = initVec3(0.f,0.f,0.f);
            for( m = 0; m < numNeighbours; m++ )
            {
                vec3 tmp;
                if( m+1 == numNeighbours )
                    tmp = cross( offset[m],offset[0]);
                else
                    tmp = cross( offset[m],offset[m+1]);
                sum.x += tmp.x;
                sum.y += tmp.y;
                sum.z += tmp.z;
            }
            vec3 result = normalize( sum );
           const int index = i*width+j;
           normalMap[index].x = result.x; normalMap[index].y = result.y; normalMap[index].z = result.z;
     }
 }

 /**
  * Reduce the overshooting phenomenon when the wave enters a shallow region
  */
 __global__ void overshootingReduction( const float* depthMap, float* nextDepthMap, const float* heightMap,
                                        const float dx, const int width, const int height )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 0 && i < height && j >= 0 && j < width )
     {
         float hij = map2Dread( depthMap, i,j,width );
         // Exclude the border
         if( i >= 1 && i < height-1 && j >= 1 && j < width-1 )
         {
             // 2.2 section
             const float alpha = 0.3;
             // n(i,j)
             float n = map2Dread( heightMap, i,j,width );
             // n(i-1,j)
             float n2 = map2Dread( heightMap, i-1,j,width );
             // n(i+1,j)
             float n3 = map2Dread( heightMap, i+1,j,width );
             // n(i,j-1)
             float n4 = map2Dread( heightMap, i, j-1, width );
             // n(i,j+1)
             float n5 = map2Dread( heightMap, i, j+1, width );
             float value;
             float nextD = hij;
             float lamda = 2*dx;
             if( n - n2 >lamda && n > n3  )
             {
                 value = alpha*( cudaMax( 0.f, 0.5*( hij + map2Dread(depthMap,i+1,j,width) ) ) - hij );
                 nextD += value;
             }
             if( n - n3 > lamda && n > n2 )
             {
                 value = alpha*( cudaMax( 0.f, 0.5*( hij + map2Dread(depthMap,i-1,j,width) ) ) - hij );
                 nextD += value;
             }
             if( n - n4 > lamda && n > n5 )
             {
                 value = alpha*( cudaMax( 0.f, 0.5*( hij + map2Dread(depthMap, i,j+1,width ) ) ) - hij );
                 nextD += value;
             }
             if( n - n5 > lamda && n > n4 )
             {
                 value = alpha*( cudaMax( 0.f, 0.5*( hij + map2Dread(depthMap, i, j- 1, width) ) ) - hij );
                 nextD += value;
             }
             map2Dwrite( nextDepthMap, i,j, nextD,width );
         }
         else
         {
             // Just copy
             map2Dwrite( nextDepthMap, i,j, hij,width );
         }
     }
 }

/**
 * Initialize the particle positions field
 **/
__global__ void initParticlePositionsCUDA( vec3* positionsMap, int numParticles )
{
    //index of current vector
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= 0 && i < numParticles ){
        positionsMap[i].x = 0;
        positionsMap[i].y = -1;
        positionsMap[i].z = 0;
    }
}

/**
 * Initialize the particle velocities field
 **/
__global__ void initParticleVelocitiesCUDA( vec3* velocitiesMap, int numParticles )
{
    //index of current vector
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= 0 && i < numParticles ){
        velocitiesMap[i].x = 0;
        velocitiesMap[i].y = 0;
        velocitiesMap[i].z = 0;
    }
}

/**
 * Update the particle positions and velocities fields
 **/
__global__ void updateParticleValuesCUDA( vec3* positionsMap, vec3* velocitiesMap, float accX, float accY, float accZ, float dt, int numParticles )
{
    //index of current vector
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= 0 && i < numParticles ){
        //check if particle is active
        if(positionsMap[i].y >= 0){
            //update position vector
            positionsMap[i].x = positionsMap[i].x +
                    (velocitiesMap[i].x * dt) +
                    (accX * dt * dt);
            positionsMap[i].y = positionsMap[i].y +
                    (velocitiesMap[i].y * dt) +
                    (accY * dt * dt);
            positionsMap[i].z = positionsMap[i].z +
                    (velocitiesMap[i].z * dt) +
                    (accZ * dt * dt);

            //update velocity vector
            velocitiesMap[i].x = velocitiesMap[i].x + (accX * dt);
            velocitiesMap[i].y = velocitiesMap[i].y + (accY * dt);
            velocitiesMap[i].z = velocitiesMap[i].z + (accZ * dt);
        }
    }
}

// TODO: CHECK THIS
// TODO: CONCERNS:
//  1. no concurrency on depth map, might lose contribution
/**
 * Intersect the particles with the height and depth fields
 **/
__global__ void updateParticleValuesCUDA( vec3* positionsMap, float* heightMap, float* depthMap,
                                          const int width, const int height,
                                          const float halfDomain, const float mdxInv, float heightChange, int numParticles )
{
    //index of current vector
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if( i >= 0 && i < numParticles ){
        //check if particle is active
        if(positionsMap[i].y >= 0){
            //find grid positions x and z
            float lenX = (positionsMap[i].x + halfDomain) * mdxInv;
            float lenZ = (positionsMap[i].z + halfDomain) * mdxInv;
            int x = (int) cudaMin(width - 1, cudaMax(0.0, round(lenX)));
            int z = (int) cudaMin(height - 1, cudaMax(0.0, round(lenZ)));

            //check if position y < heightMap
            float eta = map2Dread( heightMap, z, x, width );
            if(eta >= positionsMap[i].y){
                float hxz = map2Dread( depthMap, z, x, width );
                map2Dwrite( depthMap, z, x, hxz + heightChange, width );

                positionsMap[i].y = -1;
            }
        }
    }
}

// TODO: CHECK THIS
/**
 *  clamp the depth (min 0)
 */
__global__ void clampDepthCUDA( float* depthMap, const int width, const int height )
{
    int i = blockDim.y*blockIdx.y +threadIdx.y;
    int j = blockDim.x*blockIdx.x +threadIdx.x;

    if( i >= 0 && i < height && j >= 0 && j < width )
    {
        float depth = cudaMax(0.f, map2Dread( depthMap, i, j, width ));
        map2Dwrite( depthMap, i, j, depth, width );
    }
}

// TODO: CHECK THIS
/**
 * clamp the velocity (max velocity clamp)
 **/
__global__ void clampFieldCUDA( float* deviceMap, float velocityClamp, int width, int height )
{
    int i = blockDim.y*blockIdx.y + threadIdx.y;
    int j = blockDim.x*blockIdx.x + threadIdx.x;
    if( i >= 0&& i < height&&j >= 0 && j < width )
    {
        float value = cudaMin(velocityClamp, map2Dread( deviceMap, i, j, width ));
        map2Dwrite( deviceMap, i, j, value, width );
    }
}

//TODO: CHECK THIS
/**
 *  Initialize the sigma and gamma fields
 */
__global__ void initSigmaGammaCUDA( float* sigmaMap, float* gammaMap, const float dampeningRegion,
                                    const float quadraticA, const float quadraticB, const float quadraticC,
                                    const int width, const int height )
{
    int i = blockDim.y*blockIdx.y +threadIdx.y;
    int j = blockDim.x*blockIdx.x +threadIdx.x;
    if( i >= 0 && i < height && j >= 0 && j < width )
    {
        if(i < dampeningRegion || i >= height - dampeningRegion ||
                j < dampeningRegion || j >= width - dampeningRegion){
            //compute horizontal and vertical distances
            float iDistance = 0;
            float jDistance = 0;
            if(i < dampeningRegion){
                iDistance = dampeningRegion - i;
            } else if(i >= height - dampeningRegion){
                iDistance = i - (height - dampeningRegion - 1);
            }

            if(j < dampeningRegion){
                jDistance = dampeningRegion - j;
            } else if(j >= width - dampeningRegion){
                jDistance = j - (width - dampeningRegion - 1);
            }

            iDistance /= (float)dampeningRegion;
            jDistance /= (float)dampeningRegion;

            //distance
            float distance = sqrt((iDistance * iDistance) + (jDistance * jDistance));

            //quadratic function
            float value = (quadraticA * distance * distance) + (quadraticB * distance) + quadraticC;

            //initialize to value
            map2Dwrite( sigmaMap, i, j, value, width );
            map2Dwrite( gammaMap, i, j, value, width );
        } else {
            //initialize to 0
            map2Dwrite( sigmaMap, i, j, 0.0f, width );
            map2Dwrite( gammaMap, i, j, 0.0f, width );
        }
    }
}

/**
 *  Initialize the phi and psi fields
 */
//TODO: CHECK THIS
__global__ void initPhiPsiCUDA( float* phiMap, float* psiMap, const int width, const int height )
{
    int i = blockDim.y*blockIdx.y +threadIdx.y;
    int j = blockDim.x*blockIdx.x +threadIdx.x;
    if( i >= 0 && i < height && j >= 0 && j < width )
    {
        //initialize to 0
        map2Dwrite( phiMap, i, j, 0.0f, width );
        map2Dwrite( psiMap, i, j, 0.0f, width );
    }
}

/**
 *  Dampen the waves (dampen the depth field), update wave dampening data structures
 */
//TODO: CHECK THIS
__global__ void dampenWavesCUDA( float* depthMap, float* heightMap, float* velUMap, float* velWMap,
                                 float* sigmaMap, float* gammaMap, float* phiMap, float* psiMap,
                                 float dampeningRegion, float hRest, float dt, float dxInv, float lambdaUpdate, float lambdaDecay,
                                 const int width, const int height, const int uwidth, const int uheight, const int wwidth, const int wheight )
{
    int i = blockDim.y*blockIdx.y +threadIdx.y;
    int j = blockDim.x*blockIdx.x +threadIdx.x;

    //TODO: CORRECT ITERATION AREA?
    if( i >= 1 && i < height - 1 && j >= 1 && j < width - 1 )
    {
        if(i < dampeningRegion || i >= height - dampeningRegion ||
                j < dampeningRegion || j >= width - dampeningRegion){
            // current values
            float currH = map2Dread( heightMap, i, j, width );
            float currDepth = map2Dread( depthMap, i, j, width );

            float currVelU = map2Dread( velUMap, i, j, uwidth );
            float currVelUDec = map2Dread( velUMap, i - 1, j, uwidth );
            float currVelW = map2Dread( velWMap, i, j, wwidth );
            float currVelWDec = map2Dread( velWMap, i, j - 1, wwidth );

            float currSigma = map2Dread( sigmaMap, i, j, width );
            float currSigmaInc = map2Dread( sigmaMap, i + 1, j, width );
            float currGamma = map2Dread( gammaMap, i, j, width );
            float currGammaInc = map2Dread( gammaMap, i, j + 1, width );
            float currPhi = map2Dread( phiMap, i, j, width );
            float currPsi = map2Dread( psiMap, i, j, width );

            // Equation 10
            // h(i,j) += ((-sigma(i,j) * (h(i,j) - hRest)) + phi(i,j)) * delta_t
            // Equation 21
            // h(i,j) += ((-gamma(i,j) * (h(i,j) - hRest)) + psi(i,j)) * delta_t
            float eq10 = ((-currSigma * (currH - hRest)) + currPhi) * dt;
            float eq21 = ((-currGamma * (currH - hRest)) + currPsi) * dt;
            map2Dwrite( depthMap, i, j, currDepth + eq10 + eq21, width );

            // Equation 11
            // u(i+0.5,j) += -0.5 * (sigma(i+1,j) + sigma(i,j)) * u(i+0.5,j) * delta_t
            float eq11 = -0.5 * (currSigmaInc + currSigma) * currVelU * dt;
            map2Dwrite( velUMap, i, j, currVelU + eq11, uwidth);

            // Equation 22
            // w(i,j+0.5) += -0.5 * (gamma(i,j+1) + gamma(i,j)) * w(i,j+0.5) * delta_t
            float eq22 = -0.5 * (currGammaInc + currGamma) * currVelW * dt;
            map2Dwrite( velWMap, i, j, currVelW + eq22, wwidth);

            //TODO: reread u and w????

            // Equation 12
            // phi(i,j) += -LAMBDA_UPDATE * sigma(i,j) * ((w(i,j+0.5) - w(i,j-0.5)) / delta_x) * delta_t
            // Equation 13
            // phi(i,j) *= LAMBDA_DECAY
            float eq12 = -lambdaUpdate * currSigma * (currVelW - currVelWDec) * dxInv * dt;
            map2Dwrite( phiMap, i, j, (currPhi + eq12) * lambdaDecay, width );

            // Equation 23
            // psi(i,j) += -LAMBDA_UPDATE * gamma(i,j) * ((u(i+0.5,j) - u(i-0.5,j)) / delta_x) * delta_t
            // Equation 24
            // psi(i,j) *= LAMBDA_DECAY
            float eq23 = -lambdaUpdate * currGamma * (currVelU - currVelUDec) * dxInv * dt;
            map2Dwrite( psiMap, i, j, (currPsi + eq23) * lambdaDecay, width );
        }
    }
}

 // Review passed
/**
 * @brief initGrid Initialize our grid
 * @param girdSize The gridSize
 * @param terrainMap The terrainMap from host
 */
void initGridGPU( const int hostGridSize, const float hostdx, const float halfdm, const float* hostTerrainMap )
{
    gridSize = hostGridSize;
    halfDomain = halfdm;
    uwidth = gridSize + 1;
    uheight = gridSize;
    wwidth = gridSize;
    wheight = gridSize+1;
    mapdx = hostdx;
    mapdxInv = 1.f/mapdx;

    int width = gridSize;
    int height = gridSize;
    // Firstly backup the terrain's heightMap
    error = cudaMalloc(&deviceTerrainMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( deviceTerrainMap );
    error = cudaMemcpy( deviceTerrainMap, hostTerrainMap, width*height*sizeof(float), cudaMemcpyHostToDevice );
    checkCudaError( error );

    // Malloc heightMap
    error = cudaMalloc(&deviceHeightMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( deviceHeightMap );

    // Malloc depthMap
    error = cudaMalloc(&deviceDepthMap,  width*height*sizeof(float));
    checkCudaError( error );
    check1DNotNull( deviceDepthMap );

    // Malloc nextDepthMap
    error = cudaMalloc(&deviceNextDepthMap,  width*height*sizeof(float));
    checkCudaError( error );
    check1DNotNull( deviceNextDepthMap );

    // Malloc the normapMap
    error = cudaMalloc(&deviceNormalMap,  width*height*sizeof(vec3));
    checkCudaError( error );
    check1DNotNull( deviceNormalMap );

    // Malloc the paintMap
    error = cudaMalloc(&devicePaintMap,  width*height*sizeof(vec3));
    checkCudaError( error );
    check1DNotNull( devicePaintMap );

    width = uwidth;
    height = uheight;
    // Malloc velocityUMap
    error = cudaMalloc(&deviceVelocityUMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( deviceVelocityUMap );
    // Malloc nextVelocityUMap
    error = cudaMalloc(&deviceNextVelocityUMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( deviceNextVelocityUMap );

    width  =wwidth;
    height = wheight;
    // Malloc velocityWMap
    error = cudaMalloc(&deviceVelocityWMap, width*height*sizeof(float));
    checkCudaError( error );
    check1DNotNull( deviceVelocityWMap );
    // Malloc velocityWMap
    error = cudaMalloc(&deviceNextVelocityWMap, width*height*sizeof(float));
    checkCudaError( error );
    check1DNotNull( deviceNextVelocityWMap );

    // initialize the depth map
    dim3 threadsPerBlock(blockSizeX,blockSizeY);
    int blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    int blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    dim3 blocksPerGrid(blockPerGridX,blockPerGridY);
    initDepthCUDA<<<blocksPerGrid,threadsPerBlock>>>(
                                                   deviceDepthMap,deviceTerrainMap,gridSize,gridSize
                                                   );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
    // Copy the depth field to initialize the next depth map
    cudaMemcpy(deviceNextDepthMap,deviceDepthMap,gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToDevice );

    // Initialize the normal map
    initNormalCUDA<<<blocksPerGrid,threadsPerBlock>>>(
                                                   deviceNormalMap,gridSize,gridSize
                                                   );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    // Initialize velocity U map
    blockPerGridX = (uwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (uheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    initFieldCUDA<<<blocksPerGrid, threadsPerBlock>>>(deviceVelocityUMap, uwidth, uheight );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
    // Copy the velocity U map to initialize next velocity U map
    cudaMemcpy(deviceNextVelocityUMap,deviceVelocityUMap,(uwidth)*uheight*sizeof(float), cudaMemcpyDeviceToDevice );


    //checkInitializedDeviceField( deviceVelocityUMap, gridSize+1, gridSize );
    //checkInitializedDeviceField( deviceNextVelocityUMap, gridSize+1, gridSize );


    // Initialize velocityW
    blockPerGridX = (wwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (wheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    initFieldCUDA<<<blocksPerGrid, threadsPerBlock>>>(deviceVelocityWMap, wwidth, wheight );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
    // Copy the velocity W map to initialize next velocity W map
    cudaMemcpy(deviceNextVelocityWMap,deviceVelocityWMap,(wwidth)*wheight*sizeof(float), cudaMemcpyDeviceToDevice );


    //checkInitializedDeviceField( deviceVelocityWMap, wwidth, wheight );
    //checkInitializedDeviceField( deviceNextVelocityWMap, wwidth, wheight );


    updateHeightCUDA<<<blocksPerGrid,threadsPerBlock>>>(
                                                   deviceHeightMap,deviceDepthMap,deviceTerrainMap,
                                                   gridSize,gridSize
                                                   );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    updatePaintCUDA<<<blocksPerGrid,threadsPerBlock>>>( devicePaintMap, deviceHeightMap, halfDomain, mapdx, gridSize,gridSize);
    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

// Review passed
/**
 * @brief updateFluidGPU Update function interface
 */
void updateFluidGPU( const float dt )
{
    /**
     * Advect the depth
     */
    dim3 threadsPerBlock(blockSizeX,blockSizeY);
    int blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    int blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    dim3 blocksPerGrid(blockPerGridX,blockPerGridY);
    advectDepthCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceDepthMap, deviceNextDepthMap, deviceVelocityUMap,deviceVelocityWMap,
                     gridSize,gridSize,dt,mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
    /**
     * Copy back the depth
     */
    cudaMemcpy( deviceDepthMap, deviceNextDepthMap,
                gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToDevice);
    /**
     * Advect the velocity U
     */
    blockPerGridX = (uwidth+blockSizeX-1)/(blockSizeX);
    blockPerGridY = (uheight + blockSizeY-1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    advectVelUCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceVelocityUMap, deviceNextVelocityUMap,
                                                       deviceVelocityWMap, uwidth,uheight, dt,mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Copy back the velocity U
     */
    cudaMemcpy( deviceVelocityUMap, deviceNextVelocityUMap,
                uwidth*uheight*sizeof(float), cudaMemcpyDeviceToDevice);
    /**
     * Advect the velocity W
     */
    blockPerGridX = (wwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (wheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    advectVelWCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceVelocityWMap, deviceNextVelocityWMap,
                    deviceVelocityUMap, wwidth, wheight, dt, mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Copy back the velocity W
     */
    cudaMemcpy( deviceVelocityWMap, deviceNextVelocityWMap, wwidth*wheight*sizeof(float),cudaMemcpyDeviceToDevice );

    /**
     * Update the depth
     */
    blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    updateDepthCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceDepthMap,
                                                        deviceVelocityUMap, deviceVelocityWMap, gridSize,gridSize,dt,mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Update the height
     */
    updateHeightCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceHeightMap,
                                                         deviceDepthMap, deviceTerrainMap, gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Update the velocity U
     */
    blockPerGridX = (uwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (uheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    updateVelUCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceVelocityUMap, deviceHeightMap, deviceDepthMap, uwidth, uheight, dt, mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Update the velocity W
     */
    blockPerGridX = (wwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (wheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    updateVelWCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceVelocityWMap, deviceHeightMap, deviceDepthMap, wwidth, wheight, dt, mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Apply the boundary
     */
/*    blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    applyBoundaryCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceDepthMap, deviceHeightMap, deviceTerrainMap, gridSize,gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
*/
    overshootingReduction<<<blocksPerGrid,threadsPerBlock>>>( deviceDepthMap,
                                                              deviceNextDepthMap, deviceHeightMap, mapdx, gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
    /**
     * Copy back the buffer into deviceDepthMap
     */
    cudaMemcpy( deviceDepthMap, deviceNextDepthMap, gridSize*gridSize*sizeof(float),cudaMemcpyDeviceToDevice );

    /**
     * Apply boundary again
     **/
    applyBoundaryCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceDepthMap, deviceHeightMap, deviceTerrainMap, gridSize,gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * update the height map
     */
    updateHeightCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceHeightMap, deviceDepthMap, deviceTerrainMap, gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Compute he normal map
     */
    computeNormalCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceNormalMap, deviceHeightMap, gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Update the paint map
     */
    updatePaintCUDA<<<blocksPerGrid,threadsPerBlock>>>( devicePaintMap, deviceHeightMap, halfDomain, mapdx, gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

/**
 * @brief addDropGPU Add drop interface
 * @param posX The x coordinate
 * @param posZ The y coordinate
 * @param radius The radius
 * @param h The height added
 */
void addDropGPU(const int posX, const int posZ, const int radius, const float h )
{
    dim3 threadsPerBlock(blockSizeX,blockSizeY);
    int blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    int blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    dim3 blocksPerGrid(blockPerGridX,blockPerGridY);
    addDropCUDA<<<blocksPerGrid,threadsPerBlock>>>(
                                                   deviceDepthMap,posX,posZ,radius,h,gridSize,gridSize
                                                   );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    updateHeightCUDA<<<blocksPerGrid,threadsPerBlock>>>(
                                                   deviceHeightMap,deviceDepthMap,deviceTerrainMap,
                                                   gridSize,gridSize
                                                   );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

 /*   computeNormalCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceNormalMap, deviceHeightMap, gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);*/
}

/**
 * @brief destroyCUDAmem Destroy the cuda memory
 */
void destroyGPUmem()
{
    cudaFree( devicePaintMap );
    cudaFree( deviceNormalMap );
    cudaFree( deviceTerrainMap );
    cudaFree( deviceHeightMap );
    cudaFree( deviceDepthMap );
    cudaFree( deviceVelocityUMap );
    cudaFree( deviceVelocityWMap );
    cudaFree( deviceNextDepthMap );
    cudaFree( deviceNextVelocityUMap );
    cudaFree( deviceNextVelocityWMap );

    cudaFree( deviceParticlePositionsArray );
    cudaFree( deviceParticleVelocitiesArray );

    cudaThreadExit();
    cudaDeviceReset();
}

/**
 * @brief copyback After each update we need to copy back the map
 * @param host The target host map
 */
void copybackGPU(FieldType type, float* hostMap  )
{
    switch( type )
    {
    case HEIGHT:
    {
        error = cudaMemcpy( hostMap, deviceHeightMap, gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToHost);
        break;
    }
    case DEPTH:
    {
        error = cudaMemcpy( hostMap, deviceDepthMap, gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToHost);
        break;
    }
    case NORMAL:
    {
        error = cudaMemcpy( hostMap, deviceNormalMap, gridSize*gridSize*sizeof(vec3), cudaMemcpyDeviceToHost);
        break;
    }
    case PAINT:
    {
         error = cudaMemcpy( hostMap, devicePaintMap, gridSize*gridSize*sizeof(vec3), cudaMemcpyDeviceToHost);
         break;
    }
    case PARTICLE_POSITIONS:
    {
        error = cudaMemcpy( hostMap, deviceParticlePositionsArray, deviceNumParticles * sizeof(vec3), cudaMemcpyDeviceToHost);
        break;
    }
    case PARTICLE_VELOCITIES:
    {
        error = cudaMemcpy( hostMap, deviceParticleVelocitiesArray, deviceNumParticles * sizeof(vec3), cudaMemcpyDeviceToHost);
        break;
    }
    case SIGMA:
    {
        error = cudaMemcpy( hostMap, deviceSigmaMap, gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToHost);
        break;
    }
    case GAMMA:
    {
        error = cudaMemcpy( hostMap, deviceGammaMap, gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToHost);
        break;
    }
    case PHI:
    {
        error = cudaMemcpy( hostMap, devicePhiMap, gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToHost);
        break;
    }
    case PSI:
    {
        error = cudaMemcpy( hostMap, devicePsiMap, gridSize*gridSize*sizeof(float), cudaMemcpyDeviceToHost);
        break;
    }
    default:
    {
        assert(0);
        break;
    }
    }
    checkCudaError( error );
}

/**
 * @brief findSupportGPU Find supported CUDA device counts
 * @return True if device count is not zero
 */
bool findSupportDevice()
{
       int deviceCount = 0;

       cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

       if (error_id != cudaSuccess)
       {
           printf("cudaGetDeviceCount returned error code: %d\n-> %s\n", (int)error_id, cudaGetErrorString(error_id));
           printf("> FAILED %s sample finished, exiting...\n" );
           // I turn on the exit, it will never return
           exit(EXIT_FAILURE);
           return false;
       }
       if (deviceCount == 0)
       {
           printf("> There are no device(s) supporting CUDA\n");
           return false;
       }
       else
       {
           printf("> Found %d CUDA Capable Device(s)\n", deviceCount);
           return true;
       }
}

void initParticlesGPU(const int numParticles){
    deviceNumParticles = numParticles;

    //malloc the arrays
    error = cudaMalloc(&deviceParticlePositionsArray,  deviceNumParticles * sizeof(vec3));
    checkCudaError( error );
    check1DNotNull( deviceParticlePositionsArray );

    error = cudaMalloc(&deviceParticleVelocitiesArray,  deviceNumParticles * sizeof(vec3));
    checkCudaError( error );
    check1DNotNull( deviceParticleVelocitiesArray );

    //set up the iterator properties
    int threadsPerBlock = 256;
    int blocksPerGrid = (deviceNumParticles + threadsPerBlock - 1) / threadsPerBlock;

    //initialize positions
    initParticlePositionsCUDA<<<blocksPerGrid, threadsPerBlock>>>(
                                                    deviceParticlePositionsArray,
                                                    deviceNumParticles
                                                    );

    //initialize velocities
    initParticleVelocitiesCUDA<<<blocksPerGrid, threadsPerBlock>>>(
                                                    deviceParticleVelocitiesArray,
                                                    deviceNumParticles
                                                    );

    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

void updateParticlesGPU( const float dt, const float accX, const float accY, const float accZ ){
    //set up the iterator properties
    int threadsPerBlock = 256;
    int blocksPerGrid = (deviceNumParticles + threadsPerBlock - 1) / threadsPerBlock;

    //update positions and velocities
    updateParticleValuesCUDA<<<blocksPerGrid, threadsPerBlock>>>(
                                                    deviceParticlePositionsArray,
                                                    deviceParticleVelocitiesArray,
                                                    accX, accY, accZ, dt,
                                                    deviceNumParticles
                                                    );

    //error check
    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

// TODO: CHECK THIS
void intersectParticlesGPU( const float halfDomain, const float mdxInv, const float heightChange ){
    //set up the iterator properties
    int threadsPerBlock = 256;
    int blocksPerGrid = (deviceNumParticles + threadsPerBlock - 1) / threadsPerBlock;

    //intersect particles and velocities
    updateParticleValuesCUDA<<<blocksPerGrid, threadsPerBlock>>>(
                                                    deviceParticlePositionsArray, deviceHeightMap, deviceDepthMap,
                                                    gridSize, gridSize,
                                                    halfDomain, mdxInv, heightChange, deviceNumParticles
                                                    );

    //set up grid iterator
    dim3 threadsPerBlock2(blockSizeX,blockSizeY);
    int blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    int blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    dim3 blocksPerGrid2(blockPerGridX,blockPerGridY);

    //update height field
    updateHeightCUDA<<<blocksPerGrid2, threadsPerBlock2>>>( deviceHeightMap, deviceDepthMap, deviceTerrainMap,
                                 gridSize, gridSize );

    //error check
    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

void inputParticlesGPU( const float *particlePositions, const float *particleVelocities ){
    //copy over
    error = cudaMemcpy( deviceParticlePositionsArray, particlePositions, deviceNumParticles * sizeof(vec3), cudaMemcpyHostToDevice );
    checkCudaError( error );

    error = cudaMemcpy( deviceParticleVelocitiesArray, particleVelocities, deviceNumParticles * sizeof(vec3), cudaMemcpyHostToDevice );
    checkCudaError( error );
}

// TODO: CHECK THIS
void clampFieldsGPU( const float velocityClamp ){
    //TODO: iterate over depth field, make above 0
    dim3 threadsPerBlock(blockSizeX,blockSizeY);
    int blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    int blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    dim3 blocksPerGrid(blockPerGridX,blockPerGridY);
    //clamp depth, min value is 0
    clampDepthCUDA<<<blocksPerGrid,threadsPerBlock>>>(
                                                   deviceDepthMap, gridSize, gridSize
                                                   );
    //update height field
    updateHeightCUDA<<<blocksPerGrid, threadsPerBlock>>>( deviceHeightMap, deviceDepthMap, deviceTerrainMap,
                                 gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);


    //TODO: iterate over velocity U, make below velocity clamp
    blockPerGridX = (uwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (uheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    //clamp velocity u, max value is velocityClamp
    clampFieldCUDA<<<blocksPerGrid, threadsPerBlock>>>(deviceVelocityUMap, velocityClamp, uwidth, uheight );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    //TODO: iterate over velocity W, make below velocity clamp
    blockPerGridX = (wwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (wheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    //clamp velocity w, max value is velocityClamp
    clampFieldCUDA<<<blocksPerGrid, threadsPerBlock>>>(deviceVelocityWMap, velocityClamp, wwidth, wheight );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

//TODO: CHECK THIS
void initDampeningFieldsGPU( const int sizeDampeningRegion, const float quadraticA, const float quadraticB, const float quadraticC ){
    //set size of dampening region
    deviceDampeningRegion = sizeDampeningRegion;

    //get width and height
    int width = gridSize;
    int height = gridSize;

    //initialize arrays
    //initialize sigma
    error = cudaMalloc(&deviceSigmaMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( deviceSigmaMap );

    //initialize gamma
    error = cudaMalloc(&deviceGammaMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( deviceGammaMap );

    //initialize phi
    error = cudaMalloc(&devicePhiMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( devicePhiMap );

    //initialize psi
    error = cudaMalloc(&devicePsiMap, width*height*sizeof(float) );
    checkCudaError( error );
    check1DNotNull( devicePsiMap );

    //fill arrays
    dim3 threadsPerBlock(blockSizeX,blockSizeY);
    int blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    int blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    dim3 blocksPerGrid(blockPerGridX,blockPerGridY);

    //fill sigma and gamma
    initSigmaGammaCUDA<<<blocksPerGrid, threadsPerBlock>>>( deviceSigmaMap, deviceGammaMap,
                                                    deviceDampeningRegion, quadraticA, quadraticB, quadraticC,
                                                    gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    //fill phi and psi
    initPhiPsiCUDA<<<blocksPerGrid, threadsPerBlock>>>( devicePhiMap, devicePsiMap, gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

//TODO: CHECK THIS
void dampenWavesGPU( const float hRest, const float dt, const float dxInv, const float lambdaUpdate, const float lambdaDecay ){
    dim3 threadsPerBlock(blockSizeX,blockSizeY);
    int blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    int blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    dim3 blocksPerGrid(blockPerGridX,blockPerGridY);

    //dampen waves
    dampenWavesCUDA<<<blocksPerGrid, threadsPerBlock>>>(
                                    deviceDepthMap, deviceHeightMap, deviceVelocityUMap, deviceVelocityWMap,
                                    deviceSigmaMap, deviceGammaMap, devicePhiMap, devicePsiMap,
                                    deviceDampeningRegion, hRest, dt, dxInv, lambdaUpdate, lambdaDecay,
                                    gridSize, gridSize, uwidth, uheight, wwidth, wheight
                                    );
    //update height field
    updateHeightCUDA<<<blocksPerGrid, threadsPerBlock>>>( deviceHeightMap, deviceDepthMap, deviceTerrainMap,
                                 gridSize, gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);
}

#endif
