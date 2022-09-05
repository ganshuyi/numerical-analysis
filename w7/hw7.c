
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nr.h"
#include "nrutil.h"
#include <time.h>

int main(void)
{
	srand(time(NULL));
	long idum = -(rand() % 10);
	int i,j,nrot,k,kk,l,ll;
	float temp;
	int n = 11;
	float *d, *r, **v, **mat;

	v=matrix(1,n,1,n);
	mat=matrix(1,n,1,n);
	d=vector(1,n);
	r=vector(1,n);

	// generate symmetric matrix
	for (i=1; i<=n; i++) {
		for (j=1; j<=i; j++) {
			temp = gasdev(&idum);
			//printf("%f\n",temp);
			mat[i][j] = temp;
			mat[j][i] = temp;
		}
	}

	printf("symmetric matrix:\n");
	// print (for check)
	for (i=1; i<=n; i++) {
		printf("\n");
		for (j=1; j<=n; j++) {
			printf("%.3f\t", mat[i][j]);
		}
	}
	printf("\n");
	

	// compute eigenvalues and eigenvectors of symmetric matrix
	jacobi(mat,n,d,v,&nrot);
	
	/*
	printf("\nunsorted eigenvectors:\n\n");
	for (j=1;j<=n;j++) {
		printf("eigenvalue %d = %.6f\n",j,d[j]);
		printf("eigenvectors:\n");
		for (k=1;k<=n;k++) {
			printf("%12.6f",v[k][j]);
			if ((k % 6) == 0) printf("\n");
		}
		printf("\n\n");
	}
	*/
	
	// sort eigenvalues into descending order
	eigsrt(d,v,n);
	
	printf("\nsorted eigenvectors:\n\n");
	for (j=1;j<=n;j++) {
                printf("eigenvalue %d = %.6f\n",j,d[j]);
                printf("eigenvectors:\n");
                for (k=1;k<=n;k++) {
                        printf("%12.6f",v[k][j]);
                        if ((k % 6) == 0) printf("\n");
                }
                printf("\n\n");
        }


	free_matrix(v,1,n,1,n);
	free_matrix(mat,1,n,1,n);	
	free_vector(d,1,n);
	free_vector(r,1,n);
	return 0;
}
