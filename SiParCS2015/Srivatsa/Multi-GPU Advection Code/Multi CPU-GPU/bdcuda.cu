#include <stdio.h>
#include <geninputdata.h>
#include <timer.h>
#include <bdcuda.h>

#include "cuda_runtime.h"
#include "device_launch_parameters.h" 

__global__ void odefun_cuda(real *d_Wx, real *d_Wy, real *d_GammaHV, int *d_neighbors, int size, int stencil_size, real *d_u, real *d_v, real *d_temp_sol, real *d_step_sol)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int j, index;
	real val = 0, x_sol = 0, y_sol = 0, hv_sol = 0;
	real temp_dWx, temp_dWy, temp_dGammaHV;
	real temp_du, temp_dv;

	if (tid < size)
	{
		for (j = 0; j < stencil_size; j++)
		{
			temp_dWx = d_Wx[j * stencil_size + tid];
			temp_dWy = d_Wy[j * stencil_size + tid];
			temp_dGammaHV = d_GammaHV[j * stencil_size + tid];
			index = d_neighbors[j * stencil_size + tid] - 1;
			val = d_temp_sol[index];
			x_sol = x_sol + temp_dWx * val;
			y_sol = y_sol + temp_dWy * val;
			hv_sol = hv_sol + temp_dGammaHV * val;
		}
		
		temp_du = d_u[tid];
		temp_dv = d_v[tid];
	
		d_step_sol[tid] = (-1 * temp_du * x_sol) + (-1 * temp_dv * y_sol) + hv_sol;
	}	
}


void odefun_gpu(float time_step, real *temp_sol, real *step_sol)
{
	int i, j;

	//INITIALIZE CUDA EVENTS
	cudaEvent_t start,stop;
	float elapsedTime;

	//CREATING EVENTS
	cudaEventCreate(&start);
	cudaEventCreate(&stop);

	/* pack data for gpu */

	real *d_Wx, *d_Wy, *d_GammaHV, *d_temp_sol, *d_step_sol, *d_u, *d_v;
	real *h_Wx, *h_Wy, *h_GammaHV;
	int *h_neighbors, *d_neighbors;

	h_Wx = (real *)malloc(nodes->internal_count * stencil_size * sizeof(real));
	h_Wy = (real *)malloc(nodes->internal_count * stencil_size * sizeof(real));
	h_GammaHV = (real *)malloc(nodes->internal_count * stencil_size * sizeof(real));
	h_neighbors = (int *)malloc(nodes->internal_count * stencil_size * sizeof(int));


	for (i = 0; i < nodes->internal_count; i++)
        {
                for (j = 0; j < stencil_size; j++)
                {
                        h_Wx[nodes->internal_count * j + i] = Wx->val[i][j];
                        h_Wy[nodes->internal_count * j + i] = Wy->val[i][j];
                        h_GammaHV[nodes->internal_count * j + i] = GammaHV[i][j];
                        h_neighbors[nodes->internal_count * j + i] = nodes->neighbors[i][j];
                }
        }

	cudaMalloc(&d_Wx, nodes->internal_count * stencil_size * sizeof(real));
	cudaMalloc(&d_Wy, nodes->internal_count * stencil_size * sizeof(real));
	cudaMalloc(&d_GammaHV, nodes->internal_count * stencil_size * sizeof(real));
	cudaMalloc(&d_neighbors, nodes->internal_count * stencil_size * sizeof(int));
	cudaMalloc(&d_u, nodes->total_count * sizeof(real));
	cudaMalloc(&d_v, nodes->total_count * sizeof(real));
	cudaMalloc(&d_temp_sol, nodes->total_count * sizeof(real));
	cudaMalloc(&d_step_sol, nodes->internal_count * sizeof(real));

	cudaMemcpy(d_Wx, h_Wx, nodes->internal_count * stencil_size * sizeof(real), cudaMemcpyHostToDevice);
	cudaMemcpy(d_Wy, h_Wy, nodes->internal_count * stencil_size * sizeof(real), cudaMemcpyHostToDevice);
	cudaMemcpy(d_GammaHV, h_GammaHV, nodes->internal_count * stencil_size * sizeof(real), cudaMemcpyHostToDevice);
	cudaMemcpy(d_neighbors, h_neighbors, nodes->internal_count * stencil_size * sizeof(int), cudaMemcpyHostToDevice);
	cudaMemcpy(d_u, nodes->u, nodes->total_count * sizeof(real), cudaMemcpyHostToDevice);
	cudaMemcpy(d_v, nodes->v, nodes->total_count * sizeof(real), cudaMemcpyHostToDevice);
	cudaMemcpy(d_temp_sol, temp_sol, nodes->total_count * sizeof(real), cudaMemcpyHostToDevice);

	int threads_per_block = 256;
	cudaEventRecord(start,0);
	odefun_cuda<<<((nodes->internal_count - 1) / threads_per_block) + 1, threads_per_block>>>(d_Wx, d_Wy, d_GammaHV, d_neighbors, nodes->internal_count, stencil_size, d_u, d_v, d_temp_sol, d_step_sol);
	cudaDeviceSynchronize();
	
	//FINISH RECORDING
	cudaEventRecord(stop,0);
	cudaEventSynchronize(stop);

	//CALCULATE ELAPSED TIME
	cudaEventElapsedTime(&elapsedTime,start,stop);

	//DISPLAY COMPUTATION TIME
	//printf("Elapsed Time = %f\n",elapsedTime);
	
	cudaEventRecord(start,0);
	cudaMemcpy(step_sol, d_step_sol, nodes->internal_count * sizeof(real), cudaMemcpyDeviceToHost);
        cudaEventRecord(stop,0);
        cudaEventSynchronize(stop);
        cudaEventElapsedTime(&elapsedTime,start,stop);
        printf("Elapsed Time = %f\n",elapsedTime);
	/* unpack data from gpu */

	cudaFree(d_Wx);
	cudaFree(d_Wy);
	cudaFree(d_GammaHV);
	cudaFree(d_temp_sol);
	cudaFree(d_step_sol);
	cudaFree(d_u);
	cudaFree(d_v);

	free(h_Wx);
	free(h_Wy);
	free(h_GammaHV);
	free(h_neighbors);

	return;
}
