
#ifndef TEST_CU
#define TEST_CU

#include <cuda.h>
#include <stdio.h>
#include <assert.h>

extern "C"
void testVector();

__global__ void VecAdd(float* A, float* B, float* C, int N)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < N)
        C[i] = A[i] + B[i];
}

void testVector() {
        int N = 5;
        size_t size = N * sizeof(float);

        // Allocate input vectors h_A and h_B in host memory
        float* h_A = (float*)malloc(size);
        float* h_B = (float*)malloc(size);
        float* h_C = (float*)malloc(size);

        // Initialize input vectors.

        // Allocate vectors in device memory
        void* amp;
        float* d_A;
        cudaMalloc(&d_A, size);
        float* d_B;
        cudaMalloc(&d_B, size);
        float* d_C;
        cudaMalloc(&d_C, size);
        for( int i = 0; i < N; i++ )
        {
            h_A[i] = 2;
            h_B[i] = 3;
        }

        // Copy vectors from host memory to device memory
        cudaMemcpy(d_A, h_A, size, cudaMemcpyHostToDevice);
        cudaMemcpy(d_B, h_B, size, cudaMemcpyHostToDevice);

        // Invoke kernel
        int threadsPerBlock = 256;
//        int blocksPerGrid = (N + threadsPerBlock - 1)/(threadsPerBlock);
        int blocksPerGrid = (N+threadsPerBlock-1)/threadsPerBlock;
        VecAdd<<<blocksPerGrid, threadsPerBlock>>>(d_A, d_B, d_C, N);

        // Copy result from device memory to host memory
        // h_C contains the result in host memory
        cudaMemcpy(h_C, d_C, size, cudaMemcpyDeviceToHost);

        fflush(stdin);
        fflush(stdout);
        for( int i = 0; i < N; i++ )
        {
            printf("%f", h_C[i] );
        }
        // Free device memory
        cudaFree(d_A);
        cudaFree(d_B);
        cudaFree(d_C);
}

#endif
