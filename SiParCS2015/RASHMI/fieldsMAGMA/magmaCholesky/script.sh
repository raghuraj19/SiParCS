#!/bin/tcsh
module load R/3.0.1
module load mkl/11.2.3
module load cuda/6.5 
module load intel/15.0.3
module unload ncarcompilers/1.0

rm *.o
rm *.so

R CMD SHLIB magmaCholesky.c
R CMD SHLIB magmaCholesky_m.c
R CMD SHLIB smagmaCholesky.c
R CMD SHLIB smagmaCholesky_m.c
