#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <cuda.h>
#include <cublas.h>
#include <magma.h>
#include <magma_lapack.h>
#include <time.h>

SEXP smagmaCholeskyFinal(SEXP A, SEXP n, SEXP NB, SEXP id, SEXP zeroTri, SEXP lowerTri)
{
	magma_init();
//	magma_print_devices();
	
	int In, INB, ID;
	In = INTEGER_VALUE(n);
	INB = INTEGER_VALUE(NB);
	ID = INTEGER_VALUE(id);
	double *PA = NUMERIC_POINTER(A);
	float *sPA;
	clock_t t1,t2;
        double gpu_time;

	t1 = clock();
	cudaHostAlloc((void**)&sPA,In*In*sizeof(float),cudaHostAllocPortable);
	t2 = clock();
        gpu_time = (double) (t2-t1)/CLOCKS_PER_SEC; // Magma time
        printf("Memory allocation time: %f seconds \n",gpu_time);

	int i,j;
	
	printf("Copying data double to single\n");
	t1 = clock();
	for(i = 0; i < In; i++)
	{
		#pragma simd
		for(j = 0; j < In; j++)
		{
			sPA[i*In + j] =  PA[i*In + j];
		}
	}
	t2 = clock();
        gpu_time = (double) (t2-t1)/CLOCKS_PER_SEC; // Magma time
        printf("Copy time: %f seconds \n",gpu_time);
	
	magma_int_t N, status, info, max_size;
	N = In;
	status = 0;
	magma_setdevice(ID);
	
	int lTri;
	printf("Calling potrf\n");
	lTri = INTEGER_VALUE(lowerTri);
	if(lTri){
		t1=clock();
		magma_spotrf(MagmaLower, N, sPA, N, &info);
		t2 = clock();
        	gpu_time = (double) (t2-t1)/CLOCKS_PER_SEC; // Magma time
        	printf("Computation time: %f seconds \n",gpu_time);
	}
	else{
		t1=clock();
		magma_spotrf(MagmaUpper, N, sPA, N, &info);
		t2 = clock();
        	gpu_time = (double) (t2-t1)/CLOCKS_PER_SEC; // Magma time
        	printf("Computation time: %f seconds \n",gpu_time);
	}
	if(info != 0)
	{
		printf("magma_spotrf returned error %d: %s.\n", (int) info, magma_strerror(info));
	}
	
	magma_finalize();
	cublasShutdown();
	
	//caste sPA back to double and set upper or lower triangle to zero if necessary:
	int IZeroTri = INTEGER_VALUE(zeroTri);
	int zeroUTri = IZeroTri & lTri;
	int zeroLTri = IZeroTri & !lTri;
	printf("Copying data single to double\n");
	if(!IZeroTri) {
		for(i = 1; i< In; i++) {
			for(j=1; j < In; j++) {
				PA[i*In + j] = (double) sPA[i*In + j];
			}
	}
	} else if(zeroUTri) {
		for(i = 1; i< In; i++) {
	                for(j=1; j < In; j++) {
        	                if(i > j)
                	                PA[i*In + j] = 0;
				else 
                        	        PA[i*In + j] = (double) sPA[i*In + j];
                	}      
        	}
	} else {
		for(i = 1; i< In; i++) {
                        for(j=1; j < In; j++) {
                                if(i < j)
                                        PA[i*In + j] = 0;
                                else
                                        PA[i*In + j] = (double) sPA[i*In + j];
                        }
                }
	}
	
	UNPROTECT(1);
	cudaFreeHost(sPA);
	
	return(R_NilValue);
}	
