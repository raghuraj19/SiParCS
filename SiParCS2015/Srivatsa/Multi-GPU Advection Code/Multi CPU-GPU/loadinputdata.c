#include <loadinputdata.h>

void populate_neighbor_list_from_file()
{
	FILE *f;
	int i, j;

	f = fopen("neighbors.txt", "r");

	nodes->neighbors = (int **)malloc(nodes->internal_count * sizeof(nodes->neighbors));

	for (i = 0; i < nodes->internal_count; i++)
	{
		nodes->neighbors[i] = (int *)malloc(stencil_size * sizeof(int));

		for (j = 0; j < stencil_size; j++)
		{
			fscanf(f, "%d", &nodes->neighbors[i][j]);
		}
	}

	fclose(f);

	return;
}

void populate_indexes_from_file()
{
	FILE *f;
	int i;

	int *count = (int *)malloc(9 * sizeof(int));

	f = fopen("ind.txt", "r");

	fscanf(f, "%d %d %d %d %d %d %d %d %d\n", &count[0], &count[1], &count[2], &count[3], &count[4], &count[5], &count[6], &count[7], &count[8]);

	indexes = (ind *)malloc(sizeof(ind));

	indexes->num_main = count[0];

	indexes->main = (int *)malloc(count[0] * sizeof(int));
	for (i = 0; i < count[0]; i++)
		fscanf(f, "%d\n", &indexes->main[i]);

	indexes->num_top = count[1];

	indexes->top = (int *)malloc(count[1] * sizeof(int));
	for (i = 0; i < count[1]; i++)
		fscanf(f, "%d\n", &indexes->top[i]);

	indexes->num_bottom = count[2];

	indexes->bottom = (int *)malloc(count[2] * sizeof(int));
	for (i = 0; i < count[2]; i++)
		fscanf(f, "%d\n", &indexes->bottom[i]);

	indexes->num_left = count[3];

	indexes->left = (int *)malloc(count[3] * sizeof(int));
	for (i = 0; i < count[3]; i++)
		fscanf(f, "%d\n", &indexes->left[i]);

	indexes->num_right = count[4];

	indexes->right = (int *)malloc(count[4] * sizeof(int));
	for (i = 0; i < count[4]; i++)
		fscanf(f, "%d\n", &indexes->right[i]);

	indexes->num_gtop = count[5];

	indexes->g_top = (int *)malloc(count[5] * sizeof(int));
	for (i = 0; i < count[5]; i++)
		fscanf(f, "%d\n", &indexes->g_top[i]);

	indexes->num_gbottom = count[6];

	indexes->g_bottom = (int *)malloc(count[6] * sizeof(int));
	for (i = 0; i < count[6]; i++)
		fscanf(f, "%d\n", &indexes->g_bottom[i]);

	indexes->num_gleft = count[7];

	indexes->g_left = (int *)malloc(count[7] * sizeof(int));
	for (i = 0; i < count[7]; i++)
		fscanf(f, "%d\n", &indexes->g_left[i]);

	indexes->num_gright = count[8];

	indexes->g_right = (int *)malloc(count[8] * sizeof(int));
	for (i = 0; i < count[8]; i++)
		fscanf(f, "%d\n", &indexes->g_right[i]);

	fclose(f);

	return;
}

void populate_node_coordinates_from_file()
{
	FILE *f;
	int i;

	f = fopen("nodes.txt", "r");

	nodes = (nodeset *)malloc(sizeof(nodeset));

	fscanf(f, "%d %d\n", &nodes->total_count, &nodes->internal_count);

	fclose(f);

	f = fopen("x.txt", "r");

	nodes->x = (real *)malloc(nodes->total_count * sizeof(real));

	for (i = 0; i < nodes->total_count; i++)
		fscanf(f, "%lf\n", &nodes->x[i]);

	fclose(f);

	f = fopen("y.txt", "r");

	nodes->y = (real *)malloc(nodes->total_count * sizeof(real));

	for (i = 0; i < nodes->total_count; i++)
		fscanf(f, "%lf\n", &nodes->y[i]);

	fclose(f);

	return;
}

/**
void populate_diff_mats_from_file()
{
int ret_code;
int nz, rows, cols;
int i;
FILE *f;
int row, col;
real val;



if ((f = fopen("Wx.mtx", "r")) == NULL)
exit(1);

mm_read_mtx_crd_size(f, &rows, &cols, &nz);


Wx = malloc(sizeof(diff_mat));
Wx->val = (real *)malloc(nz * sizeof(real));
Wx->row = (int *)malloc(nz * sizeof(int));
Wx->col = (int *)malloc(nz * sizeof(int));

for (i = 0; i < nz; i++)
{
fscanf(f, "%d %d %lg\n", &row, &col, &val);
Wx->val[i] = val;
Wx->row[i] = row - 1;
Wx->col[i] = col - 1;
}

if (f !=stdin) fclose(f);



if ((f = fopen("Wy.mtx", "r")) == NULL)
exit(1);

mm_read_mtx_crd_size(f, &rows, &cols, &nz);


Wy = malloc(sizeof(diff_mat));
Wy->val = (real *)malloc(nz * sizeof(real));
Wy->row = (int *)malloc(nz * sizeof(int));
Wy->col = (int *)malloc(nz * sizeof(int));

for (i = 0; i < nz; i++)
{
fscanf(f, "%d %d %lg\n", &row, &col, &val);
Wy->val[i] = val;
Wy->row[i] = row - 1;
Wy->col[i] = col - 1;
}

if (f !=stdin) fclose(f);



if ((f = fopen("Whv.mtx", "r")) == NULL)
exit(1);

mm_read_mtx_crd_size(f, &rows, &cols, &nz);


Whv = malloc(sizeof(diff_mat));
Whv->val = (real *)malloc(nz * sizeof(real));
Whv->row = (int *)malloc(nz * sizeof(int));
Whv->col = (int *)malloc(nz * sizeof(int));

for (i = 0; i < nz; i++)
{
fscanf(f, "%d %d %lg\n", &row, &col, &val);
Whv->val[i] = val;
Whv->row[i] = row - 1;
Whv->col[i] = col - 1;
}

if (f !=stdin) fclose(f);

return;
}
*/

void populate_sparse_diff_mats_from_file()
{

	int i, j;
	FILE *f;

	Wx = (diff_mat *)malloc(sizeof(diff_mat));
	Wy = (diff_mat *)malloc(sizeof(diff_mat));
	Whv = (diff_mat *)malloc(sizeof(diff_mat));

	Wx->val = (real **)malloc(nodes->internal_count * sizeof(real));
	Wy->val = (real **)malloc(nodes->internal_count * sizeof(real));
	Whv->val = (real **)malloc(nodes->internal_count * sizeof(real));

	for (i = 0; i < nodes->internal_count; i++)
	{
		Wx->val[i] = (real *)malloc(stencil_size * sizeof(real));
		Wy->val[i] = (real *)malloc(stencil_size * sizeof(real));
		Whv->val[i] = (real *)malloc(stencil_size * sizeof(real));
	}

	f = fopen("Wx.txt", "r");

	for (i = 0; i < nodes->internal_count; i++)
	{
		for (j = 0; j < stencil_size; j++)
		{
			fscanf(f, "%lf", &Wx->val[i][j]);
		}
	}

	fclose(f);

	f = fopen("Wy.txt", "r");

	for (i = 0; i < nodes->internal_count; i++)
	{
		for (j = 0; j < stencil_size; j++)
		{
			fscanf(f, "%lf", &Wy->val[i][j]);
		}
	}

	fclose(f);

	f = fopen("Whv.txt", "r");

	for (i = 0; i < nodes->internal_count; i++)
	{
		for (j = 0; j < stencil_size; j++)
		{
			fscanf(f, "%lf", &Whv->val[i][j]);
		}
	}

	generate_GammaHV();

	fclose(f);

	return;
}

void populate_matlab_sol()
{
	int i;
	FILE *f;
	f = fopen("matlab_sol.txt", "r");

	matlab_sol = (real *)malloc(nodes->total_count * sizeof(real));

	for (i = 0; i < nodes->total_count; i++)
	{
		fscanf(f, "%lf", &matlab_sol[i]);
	}

	fclose(f);
	free(matlab_sol);
}
