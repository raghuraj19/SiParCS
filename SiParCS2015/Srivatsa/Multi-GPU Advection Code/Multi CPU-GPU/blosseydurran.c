#define NUM_NODES 163422
//#define NUM_NODES 42013
//#define NUM_NODES 11286
//#define NUM_NODES 3155
//#define NUM_PARTITIONS 32
//#define NUM_PARTITIONS 16
//#define NUM_PARTITIONS 4
#define NUM_PARTITIONS 8
#define STENCIL_SIZE 37
#define GHOST_FLAG_B 1
#define GHOST_FLAG_I 1000

#include <bdconfig.h>
#include <loadinputdata.h>
#include <geninputdata.h>
#include <timer.h>
#include "mpi.h"
#include "vector.h"

/* paritioning variables */

typedef struct
{
	int global_id;
	int partition_id;
	int *neighbor_list;
	int counter;
	int ghost_flag;
	int local_id;
	int ghost_partition_id;

} node;

int rank;													
int num_procs;

char *partition_file;
char *neighbor_file;
node *local_nodes;
node *global_nodes;

int local_node_count;
int internal_node_count;
int ghost_node_count;
int total_node_count;

int *internal_node_list;

vector internal_ghost_nodes_receive[NUM_PARTITIONS];
vector boundary_ghost_nodes_receive[NUM_PARTITIONS];

int **internal_ghost_nodes_send_buffer;
int **boundary_ghost_nodes_send_buffer;

int **internal_ghost_nodes_receive_buffer;
int **boundary_ghost_nodes_receive_buffer;

int boundary_ghost_partition_counters_receive[NUM_PARTITIONS];
int internal_ghost_partition_counters_receive[NUM_PARTITIONS];

int boundary_ghost_partition_counters_send[NUM_PARTITIONS];
int internal_ghost_partition_counters_send[NUM_PARTITIONS];

int boundary_update_ghost_nodes_counters_send[NUM_PARTITIONS];
int boundary_update_ghost_nodes_counters_receive[NUM_PARTITIONS];

vector boundary_update_ghost_nodes_send[NUM_PARTITIONS];
vector boundary_update_ghost_nodes_receive[NUM_PARTITIONS];

double *global_node_data;
double *local_node_data;

int *neighbor_procs_send;
int neighbor_procs_send_count;

int *neighbor_procs_receive;
int neighbor_procs_receive_count;

MPI_Request send_request[NUM_PARTITIONS], receive_request[NUM_PARTITIONS];
MPI_Status send_status[NUM_PARTITIONS], receive_status[NUM_PARTITIONS];


/* local PDE solver data */

int *local_eval_nodes;
int local_eval_node_count;
real *local_Wx;
real *local_Wy;
real *local_Whv;


/* time-stepping variables */

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

/* MPI specific functions */




/* MPI function definitions */

void parse_input_arguments(int argc, char *argv[])
{
	if (argc != 3)
	{
		if(!rank)
		{
			printf("Usage: *.exe partitions.txt neighbors.txt\n");
			fprintf(stderr, "insufficient arguments.. exiting\n");
		}
		MPI_Finalize();
		exit(1);
	}
}

void init_global_nodes()
{
	int i;

	global_nodes = (node*)malloc(NUM_NODES * sizeof(node));

	for (i = 0; i < NUM_NODES; i++)
	{
		global_nodes[i].global_id = i + 1;
		global_nodes[i].neighbor_list = (int *)malloc(STENCIL_SIZE * sizeof(int));
		global_nodes[i].counter = 0;
		global_nodes[i].ghost_flag = -1;
		global_nodes[i].ghost_partition_id = -1;
	}

	return;
}

void free_global_nodes()
{
	int i;

	for(i = 0; i < NUM_NODES; i++)
		free(global_nodes[i].neighbor_list);

	free(global_nodes);
	return;
}

void init_local_nodes()
{
	int i;
	local_nodes = (node*)malloc(local_node_count * sizeof(node));

	for (i = 0; i < local_node_count; i++)
	{
		local_nodes[i].local_id = i + 1;
		local_nodes[i].neighbor_list = (int *)malloc(STENCIL_SIZE * sizeof(int));
		local_nodes[i].counter = 0;
		local_nodes[i].ghost_flag = -1;
	}

	return;
}

void free_local_nodes()
{
	int i;

	for(i = 0; i < local_node_count; i++)
		free(local_nodes[i].neighbor_list);

	free(local_nodes);
	return;
}

void partition_domain(char *partition_file, char *neighbor_file)
{
	int i, j, k, counter;	
	populate_partition_ids(partition_file);
	populate_neighbor_lists(neighbor_file);
	update_ghost_partition_ids();
	

	internal_node_list = (int *)calloc(internal_node_count, sizeof(int));

	for (i = 0; i < NUM_PARTITIONS; i++)
		vector_init(&internal_ghost_nodes_receive[i]);

	for (i = 0; i < NUM_PARTITIONS; i++)
		vector_init(&boundary_ghost_nodes_receive[i]);

	counter = 0;

	for (j = 0; j < NUM_NODES; j++)
	{
		if (global_nodes[j].partition_id == rank)
		{
			internal_node_list[counter++] = global_nodes[i].global_id;

			for (k = 0; k < STENCIL_SIZE - 1; k++)
			{
				int neighbor_id = global_nodes[j].neighbor_list[k] - 1;
				int neighbor_pid = global_nodes[neighbor_id].partition_id;

				if (neighbor_pid != rank && global_nodes[neighbor_id].counter == 0)
				{
					global_nodes[neighbor_id].counter = 1;
					ghost_node_count++;

					if (global_nodes[neighbor_id].ghost_flag == GHOST_FLAG_I)
					{
						internal_ghost_partition_counters_receive[neighbor_pid]++;
						vector_add(&internal_ghost_nodes_receive[neighbor_pid], (void *)(global_nodes[neighbor_id].global_id));
					}
					else
					{
						boundary_ghost_partition_counters_receive[neighbor_pid]++;
						vector_add(&boundary_ghost_nodes_receive[neighbor_pid], (void *)(global_nodes[neighbor_id].global_id));
					}
				}
			}
		}
	}

	for (i = 0; i < NUM_PARTITIONS; i++)
		vector_init(&boundary_update_ghost_nodes_receive[i]);

	for(j = 0; j < NUM_NODES; j++)
	{
		int ghost_pid = global_nodes[j].ghost_partition_id;
		int pid = global_nodes[j].partition_id;

		if(ghost_pid == rank)
		{
			boundary_update_ghost_nodes_counters_receive[pid]++;
			vector_add(&boundary_update_ghost_nodes_receive[pid], (void *)(global_nodes[j].global_id));
		}
	}

	//printf("There are %d stencil centers in %d partition\n", local_eval_node_count, rank);

	// printf("The node count on partition %d is internal: %d and ghost: %d\n", rank, internal_node_count, ghost_node_count);

	/*

	for (i = 0; i < NUM_PARTITIONS; i++)
	{
		for (j = 0; j < vector_total(&internal_ghost_nodes_receive[i]); j++)
			vector_delete(&internal_ghost_nodes_receive[i], j);
	}

	for (i = 0; i < NUM_PARTITIONS; i++)
	{
		for (j = 0; j < vector_total(&boundary_ghost_nodes_receive[i]); j++)
			vector_delete(&boundary_ghost_nodes_receive[i], j);
	}

	for (i = 0; i < NUM_PARTITIONS; i++)
		vector_free(&internal_ghost_nodes_receive[i]);

	for (i = 0; i < NUM_PARTITIONS; i++)
		vector_free(&boundary_ghost_nodes_receive[i]);

	*/

	return;
}

void update_ghost_partition_ids()
{
	/* TODO: go through the list of internal ghost nodes and identify their corresponding internal nodes */
	int i;

	int index, gindex;

	/* updating top ghost nodes */

	for (i = 0; i < indexes->num_top; i++)
	{
		index = indexes->top[i];
		gindex = indexes->g_top[i];
		global_nodes[index - 1].ghost_partition_id = global_nodes[gindex - 1].partition_id;
	}

	/* updating bottom ghost nodes */

	for (i = 0; i < indexes->num_bottom; i++)
	{
		index = indexes->bottom[i];
		gindex = indexes->g_bottom[i];
		global_nodes[index - 1].ghost_partition_id = global_nodes[gindex - 1].partition_id;
	}

	/* updating left ghost nodes */

	for (i = 0; i < indexes->num_left; i++)
	{
		index = indexes->left[i];
		gindex = indexes->g_left[i];
		global_nodes[index - 1].ghost_partition_id = global_nodes[gindex - 1].partition_id;
	}


	/* updating right ghost nodes */
	for (i = 0; i < indexes->num_right; i++)
	{
		index = indexes->right[i];
		gindex = indexes->g_right[i];
		global_nodes[index - 1].ghost_partition_id = global_nodes[gindex - 1].partition_id;
	}

	return;
}

void populate_partition_ids(char *file_name)
{
	FILE *fp;
	int pid;
	int i;

	fp = fopen(file_name, "r");

	if (fp == NULL)
	{
		fprintf(stderr, "error reading file %s.. exiting \n", file_name);
		exit(1);
	}

	for (i = 0; i < NUM_NODES; i++)
	{
		if (fscanf(fp, "%d ", &pid) == EOF)
		{
			fprintf(stderr, "error reading file %s.. exiting \n", file_name);
			exit(1);
		}

		global_nodes[i].partition_id = pid;

		if(rank == pid)
			internal_node_count++;
	}

	fclose(fp);

	return;
}

void populate_neighbor_lists(char *file_name)
{
	FILE *fp;
	int i, j;

	fp = fopen(file_name, "r");

	if (fp == NULL)
	{
		fprintf(stderr, "error opening file %s.. exiting\n", file_name);
		exit(1);
	}

	for (i = 0; i < NUM_NODES; i++)
	{
		if (fscanf(fp, "%d ", &(global_nodes[i].ghost_flag)) == EOF)
		{
			fprintf(stderr, "error reading file %s.. exiting \n", file_name);
			exit(1);
		}

		for (j = 0; j < STENCIL_SIZE - 1; j++)
		{
			if (fscanf(fp, "%d ", &(global_nodes[i].neighbor_list[j])) == EOF)
			{
				fprintf(stderr, "error reading file %s.. exiting \n", file_name);
				exit(1);
			}
		}
	}

	fclose(fp);

	return;
}

void determine_neighbor_procs()
{
	int i, counter;

	neighbor_procs_receive_count = 0;

	for(i = 0; i < NUM_PARTITIONS; i++)
	{
		if(internal_ghost_partition_counters_receive[i] > 0 || boundary_ghost_partition_counters_receive[i] > 0)
		{
			neighbor_procs_receive_count++;
		}
	}

	neighbor_procs_receive = (int *)calloc(neighbor_procs_receive_count, sizeof(int));

	counter = 0;

	for(i = 0; i < NUM_PARTITIONS; i++)
	{
		if(internal_ghost_partition_counters_receive[i] > 0 || boundary_ghost_partition_counters_receive[i] > 0)
		{
			neighbor_procs_receive[counter++] = i;
		}
	}

	neighbor_procs_send_count = 0;

	for(i = 0; i < NUM_PARTITIONS; i++)
	{
		if(internal_ghost_partition_counters_send[i] > 0 || boundary_ghost_partition_counters_send[i] > 0)
		{
			neighbor_procs_send_count++;
		}
	}

	neighbor_procs_send = (int *)calloc(neighbor_procs_send_count, sizeof(int));

	counter = 0;

	for(i = 0; i < NUM_PARTITIONS; i++)
	{
		if(internal_ghost_partition_counters_send[i] > 0 || boundary_ghost_partition_counters_send[i] > 0)
		{
			neighbor_procs_send[counter++] = i;
		}
	}

	/* debugging purposes */
	/*
	if(rank == 0)
	{
		for(i = 0; i < neighbor_procs_receive_count; i++)
			printf("%d\t", neighbor_procs_receive[i]);

		printf("\n");

		for(i = 0; i < neighbor_procs_send_count; i++)
			printf("%d\t", neighbor_procs_send[i]);

		printf("\n");
	}
	*/
}

void exchange_send_counts()
{
	int *global_internal_ghost_receive = (int *)calloc(NUM_PARTITIONS * NUM_PARTITIONS, sizeof(int));
	int *global_boundary_ghost_receive = (int *)calloc(NUM_PARTITIONS * NUM_PARTITIONS, sizeof(int));

	int *global_internal_ghost_send = (int *)calloc(NUM_PARTITIONS * NUM_PARTITIONS, sizeof(int));
	int *global_boundary_ghost_send = (int *)calloc(NUM_PARTITIONS * NUM_PARTITIONS, sizeof(int));

	MPI_Gather(internal_ghost_partition_counters_receive, NUM_PARTITIONS, MPI_INT, global_internal_ghost_receive, NUM_PARTITIONS, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Gather(boundary_ghost_partition_counters_receive, NUM_PARTITIONS, MPI_INT, global_boundary_ghost_receive, NUM_PARTITIONS, MPI_INT, 0, MPI_COMM_WORLD);

	transpose_array(global_internal_ghost_receive, global_internal_ghost_send);
	transpose_array(global_boundary_ghost_receive, global_boundary_ghost_send);

	/* debugging purposes */
	/*
	if(rank == 0)
	{
		for(i = 0; i < NUM_PARTITIONS * NUM_PARTITIONS; i++)
		{
			printf("%d\t%d\t", global_internal_ghost_receive[i], global_boundary_ghost_receive[i]);
			if((i+1) % 8 == 0)
				printf("\n");
		}

		for(i = 0; i < NUM_PARTITIONS * NUM_PARTITIONS; i++)
		{
			printf("%d\t%d\t", global_internal_ghost_send[i], global_boundary_ghost_send[i]);
			if((i+1) % 8 == 0)
				printf("\n");
		}
	}
	*/

	MPI_Scatter(global_internal_ghost_send, NUM_PARTITIONS, MPI_INT, internal_ghost_partition_counters_send, NUM_PARTITIONS, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatter(global_boundary_ghost_send, NUM_PARTITIONS, MPI_INT, boundary_ghost_partition_counters_send, NUM_PARTITIONS, MPI_INT, 0, MPI_COMM_WORLD);

	//free(global_internal_ghost_receive);
	//free(global_boundary_ghost_receive);

	//free(global_internal_ghost_send);
	//free(global_boundary_ghost_send);
}

void transpose_array(int *array, int *transpose)
{
	int i, j;

	for(i = 0; i < NUM_PARTITIONS; i++)
	{
		for(j = 0; j < NUM_PARTITIONS; j++)
		{
			transpose[j * NUM_PARTITIONS + i] = array[i * NUM_PARTITIONS + j];
		}
	}

	return;
}

void build_receive_buffers()
{
	int i, j;

	internal_ghost_nodes_receive_buffer = (int **)calloc(neighbor_procs_receive_count, sizeof(*internal_ghost_nodes_receive_buffer));

	for(i = 0; i < neighbor_procs_receive_count; i++)
	{
		int pid = neighbor_procs_receive[i];
		internal_ghost_nodes_receive_buffer[i] = (int *)calloc(internal_ghost_partition_counters_receive[pid], sizeof(int));
	}	

	for(i = 0; i < neighbor_procs_receive_count; i++)
	{
		int pid = neighbor_procs_receive[i];

		for(j = 0; j < vector_total(&internal_ghost_nodes_receive[pid]); j++)
		{
			internal_ghost_nodes_receive_buffer[i][j] = (int)vector_get(&internal_ghost_nodes_receive[pid],j);
		}
	}

	return;
}

void build_send_buffers()
{
	int i;

	internal_ghost_nodes_send_buffer = (int **)calloc(neighbor_procs_send_count, sizeof(*internal_ghost_nodes_send_buffer));

	for(i = 0; i < neighbor_procs_send_count; i++)
	{
		int pid = neighbor_procs_send[i];
		internal_ghost_nodes_send_buffer[i] = (int *)calloc(internal_ghost_partition_counters_send[pid], sizeof(int));
	}

	return;
}

void exchange_receive_indexes()
{
	int i;

	for(i = 0; i < neighbor_procs_receive_count; i++)
	{
		send_request[i] = MPI_REQUEST_NULL;	
		receive_request[i] = MPI_REQUEST_NULL;	
	}

	for(i = 0; i < neighbor_procs_send_count; i++)
	{
		int pid = neighbor_procs_send[i];
      	MPI_Irecv(internal_ghost_nodes_send_buffer[i], internal_ghost_partition_counters_send[pid], MPI_INT, pid, 0,
              MPI_COMM_WORLD, &receive_request[i]);      	
    }

    for(i = 0; i < neighbor_procs_receive_count; i++)
	{
		int pid = neighbor_procs_receive[i];
		MPI_Isend(internal_ghost_nodes_receive_buffer[i], internal_ghost_partition_counters_receive[pid], MPI_INT, pid, 0,
              MPI_COMM_WORLD, &send_request[i]);		
	}

	for(i = 0; i < neighbor_procs_send_count; i++)
		MPI_Wait(&send_request[i], &send_status[i]);

	for(i = 0; i < neighbor_procs_receive_count; i++)
		MPI_Wait(&receive_request[i], &receive_status[i]);
	
	/* debugging purposes */

	/*
	if(rank == 0)
	{
		for(i = 0; i < 5 ; i++)
			printf("%d\t", internal_ghost_nodes_receive_buffer[1][i]);
		printf("\n");
	}

	if(rank == 6)
	{
		for(i = 0; i < 5 ; i++)
			printf("%d\t", internal_ghost_nodes_send_buffer[0][i]);
		printf("\n");
	}
	*/

	return;
}

/* function definitions */

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

	for (i = 1; i <= 1; i++)
	{
		for (j = 0; j < sol_len; j++)
			temp_sol[j] = actual_sol[j];	
		
		update_ghost_nodes(temp_sol);
        	generate_u(time_step);
        	generate_v(time_step);
		odefun_gpu(time_steps[i - 1], temp_sol, step_1_sol);		
		factor = (real)(time_step / 2);
		update_solution(temp_sol, actual_sol, step_1_sol, factor, sol_len);
		update_ghost_nodes(temp_sol);
        	generate_u(time_step);
	        generate_v(time_step);
		odefun_gpu(time_steps[i - 1] + factor, temp_sol, step_2_sol);
		update_solution(temp_sol, actual_sol, step_2_sol, factor, sol_len);
		update_ghost_nodes(temp_sol);
        	generate_u(time_step);
        	generate_v(time_step);	
		odefun_gpu(time_steps[i - 1] + factor, temp_sol, step_3_sol);
		factor = time_step;
		update_solution(temp_sol, actual_sol, step_3_sol, factor, sol_len);
		update_ghost_nodes(temp_sol);
       		generate_u(time_step);
        	generate_v(time_step);
		odefun_gpu(time_steps[i - 1] + factor, temp_sol, step_4_sol);
		factor = (real)(time_step / 6);
		for (j = 0; j < sol_len; j++)
			actual_sol[j] = actual_sol[j] + factor * (step_1_sol[j] + (2 * step_2_sol[j]) + (2 * step_3_sol[j]) + step_4_sol[j]);

	}

	for (j = 0; j < sol_len; j++)
	{
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

	return;
}

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

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	/* partition the domain */

	parse_input_arguments(argc, argv);

	init_global_nodes();
	partition_file = argv[1];
	neighbor_file = argv[2];

	populate_node_coordinates_from_file();
	initialize_pde_state();
	generate_time_steps();
	populate_sparse_diff_mats_from_file();
	populate_indexes_from_file();
	populate_neighbor_list_from_file();

	partition_domain(partition_file, neighbor_file);

	MPI_Barrier(MPI_COMM_WORLD);

	exchange_send_counts();

	MPI_Barrier(MPI_COMM_WORLD);

	determine_neighbor_procs();	

	MPI_Barrier(MPI_COMM_WORLD);

	build_receive_buffers();

	build_send_buffers();	

	MPI_Barrier(MPI_COMM_WORLD);
	
	start_time = get_time();
	exchange_receive_indexes();
	end_time = get_time();
	printf("Time elapsed: %llu\n", end_time - start_time);	

	MPI_Barrier(MPI_COMM_WORLD);

	/* time stepping */


	pde_solver();



	free_resources();

	fflush(stdout);
	MPI_Finalize();
	return 0;
}

void build_local_eval_nodes()
{
	/* TODO: build the local eval node indexes */

	int i;

	int counter = 0;

	local_eval_nodes = (int *)malloc(internal_node_count * sizeof(int));

	for(i = 0; i < NUM_NODES; i++)
	{
		if(global_nodes[i].partition_id == rank && global_nodes[i].ghost_flag == GHOST_FLAG_I)  
		{
			local_eval_nodes[counter++] = global_nodes[i].global_id;
		}
	}

	free(local_eval_nodes);

	return;
}


