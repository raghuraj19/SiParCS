#include <stdio.h>
#include <cuda.h>
#include <time.h> 
#include "magma.h"

int main ( int argc , char ** argv )
{
    magma_init (); // initialize Magma
    clock_t start , end;
    double gpu_time ;
    int ngpu=4;
    magma_int_t info , i, j, k;
    magma_int_t m = 8192; // a - mxm matrix
    magma_int_t mm=m*m; // size of a
    double *a,*b; // 
    int iterations=2;
    // allocate matrices on the host
     magma_dmalloc_cpu ( &a , mm ); // host memory for a
     magma_dmalloc_cpu ( &b , mm ); // host memory for a
             
    //initialze a positive definite matrix
    for(i=0; i<m; i++)
    {
            for (j=0; j<=i; j++)
            {
                     a[j*m+i] = (rand()%10) / 100.0 + 0.1;
                     a[i*m+j] = (a[j*m+i]);
		     b[i*m+j] = b[j*m+i] = a[i*m+j];
            }
	    a[i*m+i] += 100.0;//large diagonal elements 
	    b[i*m+i] = a[i*m+i];
    }                                                      
    gpu_time=0;                                                                  
    for(i=0;i<iterations;i++)
    {
	    
     	    start = clock();
            magma_dpotrf_m(ngpu,MagmaLower, m, a, m, &info);
     	    end = clock();
    	    if(info != 0)
	    {
	  	printf("Error in job!%d\n",info);
	    	goto terminate;
	    }
    	    gpu_time += ((double) end-start)/(CLOCKS_PER_SEC); // Magma  time
	    for(j=0;j<m;j++)
		for( k=0;k<=j;k++)
		{
			a[k*m+j] = a[j*m+k] = b[k*m+j];
		}
    }
    
    printf (" Elapsed Time : %g sec .\n",gpu_time*1.0/iterations );
    terminate:
    free (a); // free host memory
    magma_finalize (); // finalize Magma
    return 0;
}

