#!/bin/bash

#BSUB -n 1
#BSUB -q gpgpu
#BSUB -W 1:00
#BSUB -P NCIS0002
#BSUB -J Cholesky
#BSUB -o Cholesky%J.out
#BSUB -e Cholesky%J.err

rm *.o *.exe

export MAGMA_HOME=/home/cirrascale/magma-1.6.2

export LD_LIBRARY_PATH=${MAGMA_HOME}/lib:${LD_LIBRARY_PATH}


icpc -O3 -fPIC -DADD_ -Wall -openmp -DMAGMA_SETAFFINITY -DMAGMA_WITH_MKL -DMKL_ILP64 -Xlinker -shared -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 -c rr_dpotrf3_mgpu.cpp -I/usr/local/cuda-7.0/include -I${MAGMA_HOME}/include -I${MAGMA_HOME}/control -L/opt/intel/composer_xe_2015.3.187/mkl/lib/intel64 -L${MAGMA_HOME}/lib -lmagma -lcuda -I${MAGMA_HOME}/sparse-iter/include -I${MAGMA_HOME}/sparse-iter/control -o rr_dpotrf3_mgpu.o

icpc -O3 -fPIC -DADD_ -Wall -openmp -DMAGMA_SETAFFINITY -DMAGMA_WITH_MKL -DMKL_ILP64 -Xlinker -shared -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 -c rr_dpotrf_m.cpp -I/usr/local/cuda-7.0/include -I${MAGMA_HOME}/include -I${MAGMA_HOME}/control -L${MAGMA_HOME}/lib -L/opt/intel/composer_xe_2015.3.187/mkl/lib/intel64 -lmagma -lcuda -I${MAGMA_HOME}/sparse-iter/include -I${MAGMA_HOME}/sparse-iter/control -I${MAGMA_HOME} -o rr_dpotrf_m.o

icpc -O3 -fPIC -DADD_ -Wall -openmp -DMAGMA_SETAFFINITY -DMAGMA_WITH_MKL -DMKL_ILP64 -Xlinker -shared -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 -c  testRun.cpp -I/usr/local/cuda-7.0/include -I${MAGMA_HOME}/include -L${MAGMA_HOME}/lib -L/opt/intel/composer_xe_2015.3.187/mkl/lib/intel64 -lmagma -lcuda -I${MAGMA_HOME}/control -I${MAGMA_HOME}/sparse-iter/include -I${MAGMA_HOME}/sparse-iter/control -I${MAGMA_HOME} -o testRun.o

icpc -O3 -DADD_ -Wall -openmp -DMAGMA_SETAFFINITY -DMAGMA_WITH_MKL -DMKL_ILP64  -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 -DHAVE_CUBLAS -DMIN_CUDA_ARCH=350 rr_dpotrf_m.o  rr_dpotrf3_mgpu.o  testRun.o -o compare.exe -I/usr/local/cuda-7.0/include -I${MAGMA_HOME}/include -L${MAGMA_HOME}/lib  -L/opt/intel/composer_xe_2015.3.187/mkl/lib/intel64 -L/usr/local/cuda-7.0/lib64 -lcudart  -lcuda -lmkl_intel_ilp64 -lmagma  -I${MAGMA_HOME}/control -I${MAGMA_HOME}/sparse-iter/include -I${MAGMA_HOME}/sparse-iter/control -I${MAGMA_HOME}


./compare.exe


