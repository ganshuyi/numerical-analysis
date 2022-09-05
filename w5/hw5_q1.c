
#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define NP 20
#define MP 20
#define MAXSTR 80

int main(void)
{
	int j,k,l,m,n,*indx;
	float wmax,wmin,*w,*x,*c, p, *xn, *y;
	float **a,**b,**d,**u,**v,**inv;
	FILE *fp;

	indx=ivector(1,NP);
	w=vector(1,NP);
	x=vector(1,NP);
	y=vector(1,NP);
	xn=vector(1,NP);
	c=vector(1,NP);
	a=matrix(1,NP,1,NP);
	b=matrix(1,NP,1,MP);
	d=matrix(1,NP,1,NP);
	u=matrix(1,NP,1,NP);
	v=matrix(1,NP,1,NP);
	inv=matrix(1,NP,1,NP);

	if ((fp = fopen("lineq1.dat","r")) == NULL)
		nrerror("Data file not found\n");
	while (!feof(fp)) {
		fscanf(fp,"%d %d ",&n,&m);
		for (k=1;k<=n;k++)
			for (l=1;l<=n;l++) fscanf(fp,"%f ",&a[k][l]);
		for (k=1;k<=n;k++) fscanf(fp,"%f ",&b[k][1]);
		
		printf("Q1 initial values:\n\n");
                for (k=1;k<=n;k++) {
                        for (l=1;l<=m;l++) printf("%12.6f\t", a[k][l]);
                        printf("\n");
                }
                printf("\n");
                for (l=1;l<=m;l++) printf("%12.6f\t", b[l][1]);
                printf("\n\n");
                printf ("***********************************\n\n");

		/* copy a into u */
		for (k=1;k<=n;k++)
			for (l=1;l<=n;l++) u[k][l]=a[k][l];
		/* decompose matrix a */
		svdcmp(u,n,n,w,v);
		/* find maximum singular value */
		wmax=0.0;
		for (k=1;k<=n;k++)
			if (w[k] > wmax) wmax=w[k];
		/* define "small" */
		wmin=wmax*(1.0e-6);
		/* zero the "small" singular values */
		for (k=1;k<=n;k++)
			if (w[k] < wmin) w[k]=0.0;
		
		printf("\nUsing Singular Value Decomposition:\n\n");

		/* backsubstitute for each right-hand side vector */
		for (l=1;l<=1;l++) {
			for (k=1;k<=n;k++) c[k]=b[k][l];
			svbksb(u,w,v,n,n,c,x);
			printf("Solution vector:\n");
			for (k=1;k<=n;k++) printf("%12.6f",x[k]);
			printf("\n\n(matrix)*(solution):\n");
			for (k=1;k<=n;k++) {
				c[k]=0.0;
				for (j=1;j<=n;j++)
					c[k] += a[k][j]*x[j];
			}
			for (k=1;k<=n;k++) printf("%12.6f",c[k]);
			printf("\n\n");
		}
		printf ("***********************************\n\n");
		
		printf("Using LU Decomposition:\n\n");
		for (l=1;l<=n;l++)
			for (k=1;k<=n;k++) d[k][l]=a[k][l];

		/* Do LU decomposition */
		ludcmp(d,n,indx,&p);

		//calculate determinant
                for(j=1;j<=n;j++) p *= d[j][j];
	
		/* Solve equations for each right-hand vector */
		for (l=1;l<=n;l++) x[l]=b[l][1];
		lubksb(d,n,indx,x);

		for (k=1;k<=n;k++) {
			y[k]=b[k][1];
			xn[k]=x[k];
			for (j=1;j<=n;j++) inv[j][k]=x[k];
		}
		
		printf("Solution vector before mprove:\n");
                for (k=1;k<=n;k++){        
			printf("%12.9f",x[k]);
                        printf("\n");
                }

		for (k=1;k<=1;k++) {
                        for (l=1;l<=n;l++) x[l]=b[l][1];
                        printf("\n(matrix)*(solution):\n");
                        for (l=1;l<=n;l++) {
                                b[l][k]=0.0;
                                for (j=1;j<=n;j++)
                                        b[l][k] += (a[l][j]*xn[j]);
                        }
                        for (l=1;l<=n;l++)
                              printf("%14.9f",b[l][k]);
                }

		//mprove fx
		mprove(a,d,n,indx,y,xn);

		// print results
		printf("\n\nMatrix X after mprove():\n");
                for (k=1;k<=n;k++){
                        printf("%12.9f",xn[k]);
                        printf("\n");
                }
		for (k=1;k<=1;k++) {
                        for (l=1;l<=n;l++) x[l]=b[l][1];
                        printf("\n(matrix)*(solution):\n");
                        for (l=1;l<=n;l++) {
                                b[l][k]=0.0;
                                for (j=1;j<=n;j++)
                                        b[l][k] += (a[l][j]*xn[j]);
                        }
                        for (l=1;l<=n;l++)
                              printf("%14.9f",b[l][k]);
                        printf("\n\n*********************************\n");
                }

                printf("\nDeterminant of matrix A: %f\n", p);
                
		printf("\nInverse of matrix A does not exist.\n");
		
		printf("\nSolution cannot be calculated using Gauss-Jordan elimination because matrix is singular (det=0).\n");

		printf("\n*********************************\n");

	}
	fclose(fp);
	free_matrix(v,1,NP,1,NP);
	free_matrix(u,1,NP,1,NP);
	free_matrix(b,1,NP,1,MP);
	free_matrix(a,1,NP,1,NP);
	free_matrix(d,1,NP,1,NP);
	free_matrix(inv,1,NP,1,NP);
	free_vector(c,1,NP);
	free_vector(x,1,NP);
	free_vector(w,1,NP);
	free_vector(xn,1,NP);
	free_vector(y,1,NP);
	free_ivector(indx,1,NP);
	return 0;
}
#undef NRANSI
