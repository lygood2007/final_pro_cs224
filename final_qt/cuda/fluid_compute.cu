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

vec3* deviceNormalMap; // Normal map for GPU
float* deviceTerrainMap; // Terrain map for GPU
float* deviceHeightMap; // Height map for GPU
float* deviceDepthMap; // Depth map for GPU
float* deviceVelocityUMap; // VelocityU map for GPU
float* deviceVelocityWMap; // VelocityW map for GPU

float* deviceNextDepthMap; // Temp buffer for storing next depth map
float* deviceNextVelocityUMap; // Temp buffer for storing next velocity U map
float* deviceNextVelocityWMap; // Temp buffer for storing next velocity W map

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

// Forward decaration
void initGridGPU( const int hostGirdSize, const float* hostTerrainMap );
void addDropGPU(const int posX, const int posZ, const int radius, const float h );
void copybackGPU(float* hostHeightMap );
void destroyGPU();
void advectGPU(const float dt);
void updateFluidGPU( const float dt );


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
            map2Dwrite(depthMap, i, j, dh*dt+dep, width );
     }
 }

 // Review passed
 /**
  * Update the velocity U field
  */
 __global__ void updateVelUCUDA( float* velUMap, const float* heightMap,
                                 const int width, const int height, const float dt, const float dxInv )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >=1 && i < height - 1 && j >= 2 && j < width - 2 )
     {
         // The width of heightmap is 1 smaller than the width of velocity U
        float h1 = map2Dread( heightMap, i,j, width - 1);
        float h2 = map2Dread(  heightMap, i,j-1, width-1 );

        // Read the origin value from velUMap
        float vel = map2Dread( velUMap, i,j, width );
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
 __global__ void updateVelWCUDA( float* velWMap, const float* heightMap,
                                 const int width, const int height, const float dt, const float dxInv )
 {
     int i = blockDim.y*blockIdx.y + threadIdx.y;
     int j = blockDim.x*blockIdx.x + threadIdx.x;
     if( i >= 2 && i < height - 2 && j >= 1 && j < width - 1 )
     {
         float h1 = map2Dread( heightMap, i,j, width );
         float h2 = map2Dread(  heightMap, i-1,j, width );

         float vel = map2Dread( velWMap, i,j, width );
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

 // Review passed
/**
 * @brief initGrid Initialize our grid
 * @param girdSize The gridSize
 * @param terrainMap The terrainMap from host
 */
void initGridGPU( const int hostGridSize, const float hostdx, const float* hostTerrainMap )
{
    gridSize = hostGridSize;
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
    updateVelUCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceVelocityUMap, deviceHeightMap, uwidth, uheight, dt, mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Update the velocity W
     */
    blockPerGridX = (wwidth + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (wheight + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    updateVelWCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceVelocityWMap, deviceHeightMap, wwidth, wheight, dt, mapdxInv );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

    /**
     * Apply the boundary
     */
    blockPerGridX = (gridSize + blockSizeX - 1)/(blockSizeX);
    blockPerGridY = (gridSize + blockSizeY - 1)/(blockSizeY);
    blocksPerGrid = dim3(blockPerGridX,blockPerGridY);
    applyBoundaryCUDA<<<blocksPerGrid,threadsPerBlock>>>( deviceDepthMap, deviceHeightMap, deviceTerrainMap, gridSize,gridSize );
    error = cudaDeviceSynchronize();
    checkCudaError(error);

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
    cudaFree( deviceNormalMap );
    cudaFree( deviceTerrainMap );
    cudaFree( deviceHeightMap );
    cudaFree( deviceDepthMap );
    cudaFree( deviceVelocityUMap );
    cudaFree( deviceVelocityWMap );
    cudaFree( deviceNextDepthMap );
    cudaFree( deviceNextVelocityUMap );
    cudaFree( deviceNextVelocityWMap );
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
#endif
