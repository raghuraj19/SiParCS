#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <cuda.h>
#include <cublas.h>
#include <magma.h>
#include <magma_lapack.h>

SEXP smagmaCholeskyFinal_m(SEXP A, SEXP n, SEXP NB, SEXP zeroTri, SEXP ngpu, SEXP lowerTri)
{
	magma_init();
	int ndevices;
	ndevices = INTEGER_VALUE(ngpu);
        int idevice;
        for(idevice=0; idevice < ndevices; idevice++)
        {
                magma_setdevice(idevice);
                if(CUBLAS_STATUS_SUCCESS != cublasInit())
                {
                        printf("Error: gpu %d: cublasInit failed\n", idevice);
                        magma_finalize();
                        exit(-1);
                }
        }
//	magma_print_devices();
	
	int In, INB;
	In = INTEGER_VALUE(n);
	INB = INTEGER_VALUE(NB);
	double *PA = NUMERIC_POINTER(A);
	float *sPA = calloc(In*In, sizeof(float));
	int i,j;
	for(i = 0; i < In; i++)
        {
                for(j = 0; j < In; j++)
                {
                        sPA[i*In + j] = (float) PA[i*In + j];
                }
        }
	magma_int_t N, status, info, nGPUs;
	N = In;
	status = 0;
	nGPUs = ndevices;
	
	//INB = magma_get_dpotrf_nb(N);
//	INB = 224;
//	printf("INB = %d\n", INB);
	//ngpu = ndevices;
//	printf("ngpu = %d\n", ngpu);
	//max_size = INB*(1+N/(INB*ndevices))*INB*((N+INB-1)/INB);
//	printf("max_size = %d\n", max_size);
	//int imax_size = max_size;
	//double *dA;
	//magma_dmalloc_pinned((void**)&dA, In*In*sizeof(double));
	
	//ldda = (1+N/(INB*ndevices))*INB;
//	printf("ldda = %d\n", ldda);
	//magma_dsetmatrix_1D_row_bcyclic(N, N, PA, N, dA, ldda, ngpu, INB);
	//magma_dpotrf_mgpu(ngpu, MagmaLower, N, dA, ldda, &info);
	int lTri;
	lTri = INTEGER_VALUE(lowerTri);
	if(lTri)
		magma_spotrf_m(nGPUs, MagmaLower, N, sPA, N, &info);
	else
		magma_spotrf_m(nGPUs, MagmaUpper, N, sPA, N, &info);
	if(info != 0)
	{
		printf("magma_spotrf returned error %d: %s.\n", (int) info, magma_strerror(info));
	}
	
	//magma_dgetmatrix_1D_row_bcyclic(N, N, dA, ldda, PA, N, ngpu, INB);
	//for(dev = 0; dev < ndevices; dev++)
	//{
		//magma_setdevice(dev);
		//cudaFree(dA[dev]);
	//}
	magma_finalize();
	cublasShutdown();
	
	//caste sPA back to double and set upper or lower triangle to zero if necessary:
	int IZeroTri = INTEGER_VALUE(zeroTri);
        int zeroUTri = IZeroTri & lTri;
        int zeroLTri = IZeroTri & !lTri;
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
	free(sPA);
	return(R_NilValue);
}
	
	
		
	
