//#define My_PI 3.141592653589793

#include <bdconfig.h>
#include <loadinputdata.h>
#include <geninputdata.h>
#include <timer.h>

real h = 0.005;                                                                                                 /* approximate distance between nodes */
int total_time = 1;                                                                                             /* total time length of the simulation */
int stencil_size = 37;
real center_x1 = 0.5;
real center_y1 = 0.5;
real center_x2 = 0.3;
real center_y2 = 0.5;
int K = 3;                                                                                                      /* order of the Laplacian operator for the HV, i.e. Laplacian^K */
real gamma_rbf = (real)1 / (1 << 9); 
real time_step = 0;
int num_time_steps = 0;
real GammaHV_rbf = 0;
long long start_time, end_time;

/* function declarations */


void initialize_pde_state();
void free_resources();
void update_ghost_nodes();
void update_solution(real *temp_sol, real *actual_sol, real *step_sol, real factor, int sol_len);
void pde_solver();
void odefun(float time_step, real *temp_sol, real *step_sol);
//void odefun_gpu(float time_step, real *temp_sol, real *step_sol);

/*
__global__ void odefun_cuda(real *d_Wx, real *d_Wy, real *d_GammaHV, int *d_neighbors, int size, int stencil_size, real *d_u, real *d_v, real *d_temp_sol, real *d_step_sol)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int j, index;
	real val, x_sol, y_sol, hv_sol;
	x_sol = 0;
	y_sol = 0;
	hv_sol = 0;

	if (tid < size)
	{
		for (j = 0; j < stencil_size; j++)
		{
			index = d_neighbors[tid * stencil_size + j] - 1;
			val = d_temp_sol[index];
			x_sol = x_sol + d_Wx[tid * stencil_size + j] * val;
			y_sol = y_sol + d_Wy[tid * stencil_size + j] * val;
			hv_sol = hv_sol + d_GammaHV[tid * stencil_size + j] * val;
		}

		d_step_sol[tid] = (-1 * d_u[tid] * x_sol) + (-1 * d_v[tid] * y_sol) + hv_sol;
	}	
}

*/

void update_solution(real *temp_sol, real *actual_sol, real *step_sol, real factor, int sol_len)
{
	int i;

	for (i = 0; i < sol_len; i++)
		temp_sol[i] = actual_sol[i] + (factor * step_sol[i]);
}

void pde_solver()
{
	int i, j;

	real *temp_sol, *actual_sol;
	real *step_1_sol, *step_2_sol, *step_3_sol, *step_4_sol;
	int sol_len = nodes->total_count;
	real factor;

	actual_sol = (real *)calloc(sol_len, sizeof(real));
	temp_sol = (real *)calloc(sol_len, sizeof(real));
	step_1_sol = (real *)calloc(sol_len, sizeof(real));
	step_2_sol = (real *)calloc(sol_len, sizeof(real));
	step_3_sol = (real *)calloc(sol_len, sizeof(real));
	step_4_sol = (real *)calloc(sol_len, sizeof(real));

	x_sol = (real *)malloc(nodes->internal_count * sizeof(real));
        y_sol = (real *)malloc(nodes->internal_count * sizeof(real));
        hv_sol = (real *)malloc(nodes->internal_count * sizeof(real));


	for (i = 0; i < sol_len; i++)
		actual_sol[i] = nodes->rhopsi[i];

	for (i = 1; i <= 5; i++)
	{
		for (j = 0; j < sol_len; j++)
			temp_sol[j] = actual_sol[j];	

		//start_time = get_time();

		odefun(time_steps[i - 1], temp_sol, step_1_sol);		

		//end_time = get_time();

        	//printf("Time elapsed: %llu\n", end_time - start_time);

		factor = (real)(time_step / 2);

		update_solution(temp_sol, actual_sol, step_1_sol, factor, sol_len);

                //start_time = get_time();		

		odefun(time_steps[i - 1] + factor, temp_sol, step_2_sol);

                //end_time = get_time();

                //printf("Time elapsed: %llu\n", end_time - start_time);

		update_solution(temp_sol, actual_sol, step_2_sol, factor, sol_len);

                //start_time  = get_time();

		odefun(time_steps[i - 1] + factor, temp_sol, step_3_sol);

                //end_time = get_time();

                //printf("Time elapsed: %llu\n", end_time - start_time);

		factor = time_step;

		update_solution(temp_sol, actual_sol, step_3_sol, factor, sol_len);

                //start_time = get_time();

		odefun(time_steps[i - 1] + factor, temp_sol, step_4_sol);

                //end_time = get_time();

		//printf("Time elapsed: %llu\n", end_time - start_time);

		factor = (real)(time_step / 6);

		for (j = 0; j < sol_len; j++)
			actual_sol[j] = actual_sol[j] + factor * (step_1_sol[j] + (2 * step_2_sol[j]) + (2 * step_3_sol[j]) + step_4_sol[j]);

	}

	for (j = 0; j < sol_len; j++)
	{
		//printf("%.10lf\n", actual_sol[j]);
		nodes->rhopsi[j] = actual_sol[j];
	}

	free(x_sol);
	free(y_sol);
	free(hv_sol);

	free(temp_sol);
	free(actual_sol);
	free(step_1_sol);
	free(step_2_sol);
	free(step_3_sol);
	free(step_4_sol);

	return;
}

void update_ghost_nodes(real *temp_sol)
{
	int i;
	int index, gindex;

	/* updating top ghost nodes */

	for (i = 0; i < indexes->num_top; i++)
	{
		index = indexes->top[i];
		gindex = indexes->g_top[i];
		temp_sol[gindex - 1] = temp_sol[index - 1];
	}

	/* updating bottom ghost nodes */

	for (i = 0; i < indexes->num_bottom; i++)
	{
		index = indexes->bottom[i];
		gindex = indexes->g_bottom[i];
		temp_sol[gindex - 1] = temp_sol[index - 1];
	}

	/* updating left ghost nodes */

	for (i = 0; i < indexes->num_left; i++)
	{
		index = indexes->left[i];
		gindex = indexes->g_left[i];
		temp_sol[gindex - 1] = temp_sol[index - 1];
	}


	/* updating right ghost nodes */
	for (i = 0; i < indexes->num_right; i++)
	{
		index = indexes->right[i];
		gindex = indexes->g_right[i];
		temp_sol[gindex - 1] = temp_sol[index - 1];
	}

	return;
}

void odefun(float time_step, real *temp_sol, real *step_sol)
{
	int i, j, val;

	update_ghost_nodes(temp_sol);

	generate_u(time_step);
	generate_v(time_step);

	start_time = get_time();
	for (i = 0; i < nodes->internal_count; i++)
	{
		x_sol[i] = 0;
		y_sol[i] = 0;
		hv_sol[i] = 0;
		int index;
		for (j = 0; j < stencil_size; j++)
		{
			index = nodes->neighbors[i][j] - 1;
			val = temp_sol[index];
			x_sol[i] += Wx->val[i][j] * val;
			y_sol[i] += Wy->val[i][j] * val;
			hv_sol[i] += GammaHV[i][j] * val;
		}

		step_sol[i] = ((-1) * nodes->u[i] * x_sol[i]) + ((-1) * nodes->v[i] * y_sol[i]) + hv_sol[i];
	}
	end_time = get_time();
	printf("Time elapsed: %llu\n", end_time - start_time);

	return;
}

/*

void odefun_gpu(float time_step, real *temp_sol, real *step_sol)
{
	int i, j;

	update_ghost_nodes(temp_sol);

	generate_u(time_step);
	generate_v(time_step);

	

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
			h_Wx[i * stencil_size + j] = Wx->val[i][j];
			h_Wy[i * stencil_size + j] = Wy->val[i][j];
			h_GammaHV[i * stencil_size + j] = GammaHV[i][j];
			h_neighbors[i * stencil_size + j] = nodes->neighbors[i][j];
		}
	}

        //gettimeofday(&tval_after, NULL);

        //timersub(&tval_after, &tval_before, &tval_result);

        //printf("Time elapsed: %llu.%06ld\n", (long int)tval_result.tv_sec, (long int)tval_result.tv_usec);

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

	odefun_cuda<<<((nodes->internal_count - 1) / threads_per_block) + 1, threads_per_block>>>(d_Wx, d_Wy, d_GammaHV, d_neighbors, nodes->internal_count, stencil_size, d_u, d_v, d_temp_sol, d_step_sol);
	cudaDeviceSynchronize();
	cudaMemcpy(step_sol, d_step_sol, nodes->internal_count * sizeof(real), cudaMemcpyDeviceToHost);

	//unpack data from gpu

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

*/

void initialize_pde_state()
{
	time_step = h / 20;
	num_time_steps = total_time / time_step;
	nodes->r = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->th = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->uth = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->u = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->v = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->rtilde = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->psi = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->rho = (real *)malloc(nodes->total_count * sizeof(real));
	nodes->rhopsi = (real *)malloc(nodes->total_count * sizeof(real));

	generate_r();
	generate_th();
	generate_rtilde();
	generate_rho();
	generate_psi();

	int *ii = (int *)malloc(nodes->total_count * sizeof(int));
	int ii_count = get_rtilde_indexes(ii);
	update_psi(ii, ii_count);

	generate_rhopsi();

	free(ii);

	/* Debug */

	//printf("%d\n", num_time_steps);
	//printf("%f\n", time_step);

	return;
}

void free_resources()
{

	int i;

	free(time_steps);

	for (i = 0; i < nodes->internal_count; i++)
		free(GammaHV[i]);
	free(GammaHV);

	for (i = 0; i < nodes->internal_count; i++)
		free(Wx->val[i]);
	free(Wx->val);
	free(Wx);

	for (i = 0; i < nodes->internal_count; i++)
		free(Wy->val[i]);
	free(Wy->val);
	free(Wy);

	for (i = 0; i < nodes->internal_count; i++)
		free(Whv->val[i]);
	free(Whv->val);
	free(Whv);


	free(indexes->main);
	free(indexes->top);
	free(indexes->bottom);
	free(indexes->left);
	free(indexes->right);
	free(indexes->g_top);
	free(indexes->g_bottom);
	free(indexes->g_left);
	free(indexes->g_right);
	free(indexes);


	free(nodes->r);
	free(nodes->th);
	free(nodes->uth);
	free(nodes->u);
	free(nodes->v);
	free(nodes->rtilde);
	free(nodes->psi);
	free(nodes->rho);
	free(nodes->rhopsi);
	free(nodes->x);
	free(nodes->y);
	for (i = 0; i < nodes->internal_count; i++)
		free(nodes->neighbors[i]);
	free(nodes->neighbors);
	free(nodes);

}

int main(int argc, char *argv[])
{
	int i;

	populate_node_coordinates_from_file();
	initialize_pde_state();
	generate_time_steps();
	populate_sparse_diff_mats_from_file();
	populate_indexes_from_file();
	populate_neighbor_list_from_file();

	//start_time = get_time();

	pde_solver();

	//end_time = get_time();

	//printf("Time elapsed: %llu\n", end_time - start_time);


	for (i = 0; i < nodes->total_count; i++)
	{
		//printf("%.10f\n", nodes->r[i]);
		//printf("%.10f\n", nodes->th[i]);
		//printf("%.10f\n", nodes->rtilde[i]);
		//printf("%.10f\n", nodes->psi[i]);
		//printf("%.10f\n", nodes->rhopsi[i]);
		//printf("%.9lf\n", nodes->uth[i]);
		//printf("%.10lf\n", nodes->u[i]);
	}

	free_resources();

	return 0;
}
