#include <geninputdata.h>

void generate_r()
{
	int i;

	for (i = 0; i < nodes->total_count; i++)
		nodes->r[i] = sqrt(((nodes->x[i] - center_x1) * (nodes->x[i] - center_x1)) + ((nodes->y[i] - center_y1) * (nodes->y[i] - center_y1)));

	return;
}

void generate_th()
{
	int i;

	for (i = 0; i < nodes->total_count; i++)
		nodes->th[i] = atan2(nodes->y[i] - center_y1, nodes->x[i] - center_x1);

	return;
}

void generate_rtilde()
{
	int i;

	for (i = 0; i < nodes->total_count; i++)
		nodes->rtilde[i] = 5 * sqrt(((nodes->x[i] - center_x2) * (nodes->x[i] - center_x2)) + ((nodes->y[i] - center_y2) * (nodes->y[i] - center_y2)));

	return;
}

int get_rtilde_indexes(int *ii)
{
	int i;
	int count = 0;
	for (i = 0; i < nodes->total_count; i++)
	{
		if (nodes->rtilde[i] <= 1)
		{
			ii[count] = i;
			count++;
		}
	}

	return count;
}

void generate_psi()
{
	int i;

	for (i = 0; i < nodes->total_count; i++)
		nodes->psi[i] = 0;

	return;
}

void update_psi(int *ii, int ii_count)
{
	int i, index;

	for (i = 0; i < ii_count; i++)
	{
		index = ii[i];
		nodes->psi[index] += pow((((1 + cos(My_PI * nodes->rtilde[index]))) / 2), 2);
	}

	return;
}

void generate_rho()
{
	int i;

	for (i = 0; i < nodes->total_count; i++)
		nodes->rho[i] = 1;

	return;
}

void generate_rhopsi()
{
	int i;

	for (i = 0; i < nodes->total_count; i++)
		nodes->rhopsi[i] = nodes->rho[i] * nodes->psi[i];

	return;
}

void generate_uth(real t)
{
	int i;

	//uth = @(t) 4*pi*r./T .* ( 1 - cos(2*pi*t./T) .* (1-(4*r).^6) ./ (1+(4*r).^6) );

	for (i = 0; i < nodes->total_count; i++)
	{
		nodes->uth[i] = (1 - pow((4 * nodes->r[i]), 6)) / (1 + pow((4 * nodes->r[i]), 6));
		nodes->uth[i] = nodes->uth[i] * cos((2 * My_PI * t) / total_time);
		nodes->uth[i] = 1 - nodes->uth[i];
		nodes->uth[i] = nodes->uth[i] * ((4 * My_PI * nodes->r[i]) / total_time);
	}

	return;
}

void generate_u(real t)
{
	int i;

	generate_uth(t);

	for (i = 0; i < nodes->total_count; i++)
		nodes->u[i] = nodes->uth[i] * sin(nodes->th[i]);

	return;
}

void generate_v(real t)
{
	int i;

	generate_uth(t);

	for (i = 0; i < nodes->total_count; i++)
		nodes->v[i] = -nodes->uth[i] * cos(nodes->th[i]);

	return;
}

void generate_GammaHV()
{
	int i, j;
	
	GammaHV_rbf = gamma_rbf * pow(h, (2 * K - 1));
	
	GammaHV = (real **)malloc(nodes->internal_count * sizeof(GammaHV));
	for (i = 0; i < nodes->internal_count; i++)
		GammaHV[i] = (real *)malloc(stencil_size * sizeof(real));

	for (i = 0; i < nodes->internal_count; i++)
	{
		for (j = 0; j < stencil_size; j++)
			GammaHV[i][j] = GammaHV_rbf * Whv->val[i][j];
	}

	return;
}

void generate_time_steps()
{
	int i;

	time_steps = (real *)malloc((num_time_steps + 1) * sizeof(real));

	for (i = 0; i <= num_time_steps; i++)
		time_steps[i] = i * time_step;

	/* Debug */

	for (i = 0; i <= num_time_steps; i++)
		//printf("%f\t", time_steps[i]);

		return;
}
