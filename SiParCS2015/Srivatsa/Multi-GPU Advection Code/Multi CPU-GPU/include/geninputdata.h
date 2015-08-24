#ifndef GEN_INPUT_H
#define GEN_INPUT_H

#include <bdconfig.h>

void generate_time_steps();
void generate_u(real t);
void generate_v(real t);
void generate_uth(real t);
void generate_r();
void generate_th();
void generate_rtilde();
void generate_rho();
void generate_psi();
void generate_rhopsi();
void generate_GammaHV();
int get_rtilde_indexes(int *ii);
void update_psi(int *ii, int ii_count);

#endif
