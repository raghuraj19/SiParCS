#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rdefines.h>
#include <cuda.h>
#include <cublas.h>
#include <magma.h>
#include <magma_lapack.h>

SEXP magmaCholeskyFinal(SEXP A, SEXP n, SEXP NB, SEXP id, SEXP zeroTri, SEXP lowerTri)
{
	magma_init();
//	magma_print_devices();
	
	double *h_R;
	int In, INB, ID;
	In = INTEGER_VALUE(n);
	INB = INTEGER_VALUE(NB);
	ID = INTEGER_VALUE(id);
	double *PA = NUMERIC_POINTER(A);
	int i,j;

	magma_int_t N, n2, lda, status, info, max_size;
	N=In;
   	lda = N;
   	n2 = lda*N;
  	

/*	for(i = 0; i < In; i++)
	{
		for(j = 0; j < In; j++)
		{
			printf("%.8f ", PA[i+j*In]);
		}
		printf("\n");
	}	*/

	if ( MAGMA_SUCCESS != magma_malloc_pinned( (void**) &h_R, (n2)*sizeof(double) )) {      
        fprintf( stderr, "!!!! magma_malloc_pinned failed for: %s\n", h_R ); 
        magma_finalize();                                                  
        exit(-1);   
     	}
        

        lapackf77_dlacpy( MagmaUpperLowerStr, &N, &N, PA, &lda, h_R, &lda );
	
	N = In;
	status = 0;
	magma_setdevice(ID);
	//printf("Modified by Vinay in one GPU\n");
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
		magma_dpotrf(MagmaLower, N, h_R, N, &info);
	else
		magma_dpotrf(MagmaUpper, N, h_R, N, &info);
	if(info != 0)
	{
		printf("magma_dpotrf returned error %d: %s.\n", (int) info, magma_strerror(info));
	}
		
	lapackf77_dlacpy( MagmaUpperLowerStr, &N, &N, h_R, &lda, PA, &lda );
	//magma_dgetmatrix_1D_row_bcyclic(N, N, dA, ldda, PA, N, ngpu, INB);
	//for(dev = 0; dev < ndevices; dev++)
	//{
		//magma_setdevice(dev);
		//cudaFree(dA[dev]);
	//}
	magma_free_pinned(h_R);
	magma_finalize();
	cublasShutdown();
	/*
	int IZeroTri;
        IZeroTri = INTEGER_VALUE(zeroTri);
	if(IZeroTri & lTri) {
		for(i = 1; i < In; i++)
        	{
       			for(j=0; j< i; j++)
                	{
                       		PA[i*In+j] = 0.0;
                	}
        	}
	}
	else if(IZeroTri)
		for(i = 0; i < In; i++)
                {
                        for(j=i+1; j < In; j++)
                        {
                                PA[i*In+j] = 0.0;
                        }
                }*/
	return(R_NilValue);
}
	
	
		
	
