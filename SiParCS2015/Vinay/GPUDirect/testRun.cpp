#include <common_magma.h>
#include <magma.h>

extern "C" magma_int_t
rr_dpotrf_m(
    magma_int_t ngpu,
    magma_uplo_t uplo, magma_int_t n,
    double *A, magma_int_t lda,
    magma_int_t *info);

int main()
{
	//Magma initialization
	magma_init();
	//Declaration of local variables
	double *a, *b, *dev_a, results=0;
	const int N=4098;
	int i,j;
	magma_int_t info=0, lda=N, ngpu=2;


	//Memory Allocation Segment
	magma_malloc_pinned((void**) &a,(N*N)*sizeof(double));
	magma_malloc_pinned((void**) &b,(N*N)*sizeof(double));

	
	//Generate two copies of Symmetric Positive Definite Matrix
	for(i=0;i<N;i++)
	{
		for(j=0;j<i;j++)
		{
			a[i*N+j] = 1e-9* (double)rand();
			a[j*N+i] = a[i*N+j];
			b[i*N+j] = a[i*N+j];
			b[j*N+i] = b[i*N+j];
		}
		a[i*N+i] = 1e-9*(double)rand() +  1000.0;
		b[i*N+i] = a[i*N+i];
	}

	//Call custom Magma Cholesky for obtaining results
	rr_dpotrf_m(ngpu,MagmaUpper,N,a,N,&info);

	//Call Standard Magma Cholesky for result validation
	magma_dpotrf(MagmaUpper,N,b,N,&info);
	if(info != 0)
	{
		printf("magma_dpotrf original returned error %d: %s. \n",(int) info, magma_strerror(info));
	}

	//Validate the results; Compute the RMS error value.
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			results = results + (a[i*N+j] - b[i*N+j]) * (a[i*N+j] - b[i*N+j]);
	
	//Display the results of the test
	if(results < 1e-5)
		printf("The two functions have identical results\n");
	else
		printf("The custom function had significant errors. The RMS value was %g\n",results);
	magma_free_pinned(a);
	magma_free_pinned(b);
	magma_finalize();
	return 0;

}
