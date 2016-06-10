#include "matrix.h"
#include <papi.h>

#define NUM_FLOPS 10000
#define NUM_EVENTS 3
int Events[NUM_EVENTS] = {PAPI_TOT_INS,PAPI_L1_DCM,PAPI_L2_DCM};
long_long values[NUM_EVENTS];

template <class T>
   Matrix<T> Matrix<T>::operator*(Matrix<T>& tmp){   
   long m1 = m;
   long n1 = n;

   long m2 = tmp.getM();
   long n2 = tmp.getN();

   if(n1!=m2){
      printf("Size mismatch\n");
      exit(1);
   }

   Matrix<T> C(m1,n2);

   int bSize = 64;

	if(PAPI_start_counters(Events,NUM_EVENTS) != PAPI_OK){
		printf("PAPI failed to start counters\n");
		exit(1);
	}

   clock_t time1 = clock();   
   for(int jk = 0;jk < n; jk+= bSize)
      for(int kk = 0; kk < n; kk+= bSize)
         for(int ik = 0; ik < n; ik+= bSize)
            for(int j = 0; j < bSize; j++ )
               for(int k = 0; k < bSize; k++ )
               for(int i = 0; i < bSize; i++ )
                  C.A[jk + j][ik + i] += A[jk + j][kk + k] * tmp.A[kk + k][ik + i]; 
   

   clock_t time2 = clock();   
   
	if(PAPI_read_counters(values,NUM_EVENTS) != PAPI_OK){
		printf("PAPI failed to read counters\n");
		exit(1);
	}

	printf("\n");
	printf("Instructions completed: %lld\n",values[0]);
	printf("L1 data cache misses: %lld\n",values[1]);
	printf("L2 data cache misses: %lld\n",values[2]);
	printf("\n");
   printf("Time: %f\n",(time2-time1)/(double)CLOCKS_PER_SEC);   
   
   return C;
}

