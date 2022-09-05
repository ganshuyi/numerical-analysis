
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nr.h"

#define N 100 //interval
#define NOVER2 (N/2)
#define NPTS 1000 //total data

void printHist(int *hist, int n) {
      int i, j;
      for (i = 0; i < n; i++) {
      printf("[%d] ", i);
      for ( j = 0; j < hist[i]; ++j) {
      printf("*");
      }
      printf("\n");
   }
}


int main(void)
{
	long idum = -1;
	int i,j;
	float hist[NPTS], ghist[NPTS];
	int results[NPTS/N], gresults[NPTS/N];
	
	printf("1000 samples\n\n");
	printf("uniform dist:\n");
	//generate random number using uniform distribution
	for (i=0; i<NPTS; i++) {
		hist[i] = ran1(&idum)*(2-(-3)) +(-3);
		//printf("%f\n", rn1);
	}

	//process histogram
	for (i=0; i<NPTS/N; i++) {
		for (j=0; j<NPTS;j++) {
			if ((hist[j] >= -3+i*0.5) && (hist[j] <= -2.5+i*0.5)) ++results[i];
		}
	}

	printf("\n");
	printHist(results,10);

	
	printf("gauss dev:\n");
	//generate random number using gaussian deviation
	for (i=0; i<NPTS; i++) {
                ghist[i] = 0.5 + 1.5 * gasdev(&idum);
		//printf("%f\n", rn2);
        }


	//process histogram
        for (i=0; i<10; i++) {
                for (j=0; j<NPTS;j++) {
                        if ((ghist[j] >= -NOVER2) && (ghist[j] <= NOVER2)) ++gresults[i];
                }
        }

        printf("\n");
        printHist(gresults,10);
	
	printf("\n100 samples\n\n");
        float hist100[100], ghist100[100];
        int results100[10], gresults100[10];
	
	printf("uniform dist:\n");
        //generate random number using uniform distribution
        for (i=0; i<100; i++) {
                hist100[i] = ran1(&idum)*(2-(-3)) +(-3);
                //printf("%f\n", rn1);
        }

        //process histogram
        for (i=0; i<10; i++) {
                for (j=0; j<100;j++) {
                        if ((hist100[j] >= -3+i*0.5) && (hist100[j] <= -2.5+i*0.5)) ++results100[i];
                }
        }

        printf("\n");
        printHist(results100,10);

        
        printf("gauss dev:\n");
        //generate random number using gaussian deviation
        for (i=0; i<100; i++) {
                ghist100[i] = 0.5 + 1.5 * gasdev(&idum);
                //printf("%f\n", rn2);
        }


        //process histogram
        for (i=0; i<10; i++) {
                for (j=0; j<100;j++) {
                        if ((ghist100[j] >= -NOVER2) && (ghist100[j] <= NOVER2)) ++gresults100[i];
                }
        }

        printf("\n");
        printHist(gresults100,10);


	printf("\n10000 samples\n\n");
        float hist10000[10000], ghist10000[10000];
        int results10000[10], gresults10000[10];

        printf("uniform dist:\n");
        //generate random number using uniform distribution
        for (i=0; i<10000; i++) {
                hist10000[i] = ran1(&idum)*(2-(-3)) +(-3);
                //printf("%f\n", rn1);
        }

        //process histogram
        for (i=0; i<10; i++) {
                for (j=0; j<10000;j++) {
                        if ((hist10000[j] >= -3+i*0.5) && (hist10000[j] <= -2.5+i*0.5)) ++results10000[i];
                }
        }

        printf("\n");
        printHist(results10000,10);


        printf("gauss dev:\n");
        //generate random number using gaussian deviation
        for (i=0; i<10000; i++) {
                ghist10000[i] = 0.5 + 1.5 * gasdev(&idum);
                //printf("%f\n", rn2);
        }


        //process histogram
        for (i=0; i<10; i++) {
                for (j=0; j<100;j++) {
                        if ((ghist10000[j] >= -NOVER2) && (ghist10000[j] <= NOVER2)) ++gresults10000[i];
                }
        }

        printf("\n");
        printHist(gresults10000,10);

	printf("\n100000 samples\n\n");
        float hist100000[100000], ghist100000[100000];
        int results100000[10], gresults100000[10];

        printf("uniform dist:\n");
        //generate random number using uniform distribution
        for (i=0; i<100000; i++) {
                hist100000[i] = ran1(&idum)*(2-(-3)) +(-3);
                //printf("%f\n", rn1);
        }

        //process histogram
        for (i=0; i<10; i++) {
                for (j=0; j<100000;j++) {
                        if ((hist100000[j] >= -3+i*0.5) && (hist100000[j] <= -2.5+i*0.5)) ++results100000[i];
                }
        }

        printf("\n");
        printHist(results100000,10);


        printf("gauss dev:\n");
        //generate random number using gaussian deviation
        for (i=0; i<100000; i++) {
                ghist100000[i] = 0.5 + 1.5 * gasdev(&idum);
                //printf("%f\n", rn2);
        }


        //process histogram
        for (i=0; i<10; i++) {
                for (j=0; j<100;j++) {
                        if ((ghist100000[j] >= -NOVER2) && (ghist100000[j] <= NOVER2)) ++gresults100000[i];
                }
        }

        printf("\n");
        printHist(gresults100000,10);


//	printf("%f\n",rn1);
//	printf("%f\n",rn2);
	return 0;
}
