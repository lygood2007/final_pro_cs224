
#ifndef INVERT_CU
#define INVERT_CU

#include <cuda.h>
#include <stdio.h>
#include <assert.h>

extern "C"
void invertImage(unsigned char *bits, int width, int height);

__global__ void invert(unsigned char *bits, int size)
{
    // invert one pixel
    int idx = blockIdx.x*blockDim.x + threadIdx.x;
    if(idx < size)
    {
        for(int i = 0; i < 4; i++)
        {
            bits[4*idx + i] = 0xFF - bits[4*idx + i];
        }
    }
}

void invertImage(unsigned char *bits, int width, int height) {
        unsigned char *device_bits;
        // it's a BGRA, so 4 chars per pixel
        size_t size = 4 * sizeof(unsigned char) * width * height;
        int numpixels = width*height;

        // allocate arrays on device
        cudaMalloc((void **) &device_bits, size);

        cudaMemcpy(device_bits, bits, size, cudaMemcpyHostToDevice);

        // calculation on device
        int blockSize = 32;
        int nBlocks = numpixels/blockSize + (numpixels % blockSize == 0 ? 0 : 1);
        invert <<< nBlocks, blockSize >>> (device_bits, numpixels);
        // retrieve result
        cudaMemcpy(bits, device_bits, size, cudaMemcpyDeviceToHost);
        // cleanup
        cudaFree(device_bits);
}

#endif
