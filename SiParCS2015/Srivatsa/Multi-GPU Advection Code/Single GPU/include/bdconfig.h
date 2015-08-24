#ifndef BD_CONFIG_H
#define BD_CONFIG_H

#define My_PI 3.141592653589793

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef double real;

/* state variables */

typedef struct
{
        real *x;
        real *y;
        real *u;
        real *v;
        real *uth;
        real *r;
        real *th;
        real *rtilde;
        real *rho;
        real *psi;
        real *rhopsi;
        int total_count;
        int internal_count;
        int **neighbors;
} nodeset;

typedef struct
{
        real **val;

} diff_mat;

typedef struct
{
        int *main;
        int *top;
        int *bottom;
        int *left;
        int *right;
        int *g_top;
        int *g_bottom;
        int *g_left;
        int *g_right;
        int num_main;
        int num_top;
        int num_bottom;
        int num_left;
        int num_right;
        int num_gtop;
        int num_gbottom;
        int num_gleft;
        int num_gright;

} ind;

extern real h;                                                                                                 /* approximate distance between nodes */
extern int total_time;                                                                                             /* total time length of the simulation */
real *time_steps;                                                                                               /* time step array */
extern real time_step;
extern int num_time_steps;
extern int stencil_size;
nodeset *nodes;
diff_mat *Wx, *Wy, *Whv;
ind *indexes;
extern real center_x1;
extern real center_y1;
extern real center_x2;
extern real center_y2;

extern int K;                                                                                                      /* order of the Laplacian operator for the HV, i.e. Laplacian^K */
extern real gamma_rbf;                                                                            /* tuning parameter of HV, use 2^-10 for cartesian and hex, 2^-7 for scattered */
extern real GammaHV_rbf;
real **GammaHV;

real *matlab_sol;
real *x_sol, *y_sol, *hv_sol;

extern long long start_time, end_time;



#endif
