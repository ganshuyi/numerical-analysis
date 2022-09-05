
#include <stdio.h>
#include <stdlib.h>
#include <error.h>
#define N 77

float sum (float *arr, int n) {
	float temp = 0;
	for (int i = 0; i < n; i++) {
		temp += arr[i];
	}
	return temp;
}

float sumProduct (float *arr1, float *arr2, int n) {
	float temp = 0;
	for (int i = 0; i < n; i++) {
		temp += (arr1[i] * arr2[i]);
	}
	return temp;
}

void leastSquares (float *x, float *y, float *xnew, float *ynew, int n) {
	float X, Y, Xn, Yn, X2, Y2, XXn, XYn, YXn, YYn;
       	float div, P, Q, R, S, T, U, res;
        float a, b;
	int i;

	X = sum(x, N);
        Y = sum(y, N);
        Xn = sum(xnew, N);
        Yn = sum(ynew, N);
        X2 = sumProduct(x, x, N);
        Y2 = sumProduct(y, y, N);
        XXn = sumProduct(x, xnew, N);
        XYn = sumProduct(x, ynew, N);
        YXn = sumProduct(y, xnew, N);
        YYn = sumProduct(y, ynew, N);

        div = X*X + Y*Y - N * (X2 + Y2);
	P = (Xn*X + Yn*Y - N * (XXn + YYn))/div;
        Q = (Xn*Y + Yn*X + N * (XYn - YXn))/div;
        R = (Xn - P*X - Q*Y)/N;
        S = (Yn - P*Y + Q*X)/N;
        T = (Xn*X + Yn*Y - N * (XXn + YYn))/div;
        U = (Xn*Y + Yn*X + N * (XYn - YXn))/div;

        res = 0;

        for (i=0; i<N; i++) {
		a = P*x[i] + Q*y[i] + R;
                b = T*x[i] + U*y[i] + S;
                res += ((a - xnew[i]) * (a - xnew[i]));
                res += ((b - ynew[i]) * (b - ynew[i]));
        }
	
        printf("a1: %f\ta2: %f\ta3: %f\n", P,Q,R);
        printf("a4: %f\ta5: %f\ta6: %f\n", T,U,S);
}

int main (void)
{
	int i;
	FILE *fp, *fp2, *fp3;
	float x[N], y[N], xnew[N], ynew[N];

	if ((fp = fopen("fitdata1.dat","r")) == NULL)
		perror("Data file not found\n");
	while (!feof(fp)) {
		for (i=0; i<N; i++) {
			fscanf(fp,"%f %f %f %f\n", 
					&x[i], &y[i], &xnew[i], &ynew[i]);
		}
		
		printf("File: fitdata1\n");	
		leastSquares(x, y, xnew, ynew, N);	
		printf("\n");
	}
	fclose(fp);
	
	if ((fp2 = fopen("fitdata2.dat","r")) == NULL)
                perror("Data file not found\n");
        while (!feof(fp2)) {
                for (i=0; i<N; i++) {
                        fscanf(fp2,"%f %f %f %f\n",
                                        &x[i], &y[i], &xnew[i], &ynew[i]);
                }

                printf("File: fitdata2\n");
                leastSquares(x, y, xnew, ynew, N);
                printf("\n");
        }
        fclose(fp2);
	
	if ((fp3 = fopen("fitdata3.dat","r")) == NULL)
                perror("Data file not found\n");
        while (!feof(fp3)) {
                for (i=0; i<N; i++) {
                        fscanf(fp3,"%f %f %f %f\n",
                                        &x[i], &y[i], &xnew[i], &ynew[i]);
                }

                printf("File: fitdata3\n");
                leastSquares(x, y, xnew, ynew, N);
        }
        fclose(fp3);

	return 0;
}
