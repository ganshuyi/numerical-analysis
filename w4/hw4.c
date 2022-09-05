#include <math.h> 

#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define N 1000
#define NBMAX 20

//intervals
#define X1 0.0
#define X2 400.0

#define JMAX 100 //rtbis

#define MAXIT 100 //rtflsp, rtsec

#define NMAX 100 //rtnewt

#define MAXIT2 100 //rtsafe

#define NR_END 1
#define FREE_ARG char*

// main function
float bessj0(float x)
{
	float ans;

	ans = expf(-x/200) * cos(floor(sqrt(2000 - x * x/100)/20)) - 0.01;
	return ans;
}

//zbrak
void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb)
{
	int nbb,i;
	float x,fp,fc,dx;

	nbb=0;
	dx=(x2-x1)/n;
	fp=(*fx)(x=x1);
	for (i=1;i<=n;i++) {
		fc=(*fx)(x += dx);
		if (fc*fp <= 0.0) {
			xb1[++nbb]=x-dx;
			xb2[nbb]=x;
			if(*nb == nbb) return;

		}
		fp=fc;
	}
	*nb = nbb;
}

//nrerror (nrutil.c)
void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


//vector (nrutil.c)
float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

//free_vector (nrutil.c)
void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

//bisection
float rtbis(float (*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float dx,f,fmid,xmid,rtb;

	f=(*func)(x1);
	fmid=(*func)(x2);
	if (f*fmid >= 0.0) nrerror("Root must be bracketed for bisection in rtbis");
	rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++) {
		fmid=(*func)(xmid=rtb+(dx *= 0.5));
		if (fmid <= 0.0) rtb=xmid;
		if (fabs(dx) < xacc || fmid == 0.0) {printf("iteration %d\n", j); return rtb;}
	}
	nrerror("Too many bisections in rtbis");
	return 0.0;
}

//linear intersection
float rtflsp(float (*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float fl,fh,xl,xh,swap,dx,del,f,rtf;

	fl=(*func)(x1);
	fh=(*func)(x2);
	if (fl*fh > 0.0) nrerror("Root must be bracketed in rtflsp");
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xl=x2;
		xh=x1;
		swap=fl;
		fl=fh;
		fh=swap;
	}
	dx=xh-xl;
	for (j=1;j<=MAXIT;j++) {
		rtf=xl+dx*fl/(fl-fh);
		f=(*func)(rtf);
		if (f < 0.0) {
			del=xl-rtf;
			xl=rtf;
			fl=f;
		} else {
			del=xh-rtf;
			xh=rtf;
			fh=f;
		}
		dx=xh-xl;
		if (fabs(del) < xacc || f == 0.0) 
		{printf("iteration %d\n",j); return rtf;}
	}
	nrerror("Maximum number of iterations exceeded in rtflsp");
	return 0.0;
}

//secant
float rtsec(float (*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float fl,f,dx,swap,xl,rts;

	fl=(*func)(x1);
	f=(*func)(x2);
	if (fabs(fl) < fabs(f)) {
		rts=x1;
		xl=x2;
		swap=fl;
		fl=f;
		f=swap;
	} else {
		xl=x1;
		rts=x2;
	}
	for (j=1;j<=100000;j++) {
		dx=(xl-rts)*f/(f-fl);
		xl=rts;
		fl=f;
		rts += dx;
		f=(*func)(rts);
		if (fabs(dx) < xacc || f == 0.0) 
		{printf("iteration %d\n", j);
			return rts;}
	}
	nrerror("Maximum number of iterations exceeded in rtsec");
	return 0.0;
}

//muller
float muller(float (*func)(float), float x1, float x2, float xacc)
{
	void nrerror(char error_text[]);
        int j;

	float x3, f1, f2, f3, d1, d2, h1, h2, a, b, c, x, dx, den, result;
	x3 = (x1 + x2)/2;
	
	for (j = 0;; ++j) {
		f1 = (*func)(x1);
		f2 = (*func)(x2);
		f3 = (*func)(x3);
		h1 = x2 - x1;
		h2 = x3 - x2;
		d1 = (f2 - f1)/h1;
		d2 = (f3 - f2)/h2;
		
		a = (d2 - d1)/(h2 - h1);
		b = a * h2 + d2;
		c = f3;

        	x = sqrt(b*b - 4*a*c);

		//take root closer to x2
		if (fabs(b+x) > fabs(b-x))
			den = b + x;
		else
			den = b - x;
         	dx = (-2 * c) / den;
		result = x3 + dx;
		if (fabs(dx) < xacc || j > MAXIT) 
		{	printf("iteration %d\n", j);
			return result;
		}
		x1 = x2;
         	x2 = x3;
         	x3 = result;
	}
	nrerror("Root can't be found using Muller method\n");
	return 0.0;      
}

//newton-raphson
float rtnewt(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float df,dx,f,rtn;

	rtn=0.5*(x1+x2);
	for (j=1;j<=NMAX;j++) {
		(*funcd)(rtn,&f,&df);
		dx=f/df;
		rtn -= dx;
		if ((x1-rtn)*(rtn-x2) < 0.0)
			nrerror("Jumped out of brackets in rtnewt");
		if (fabs(dx) < xacc) {printf("iteration %d\n", j); return rtn;}
	}
	nrerror("Maximum number of iterations exceeded in rtnewt");
	return 0.0;
}

//newton with bracketing
float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc)
{
	void nrerror(char error_text[]);
	int j;
	float df,dx,dxold,f,fh,fl;
	float temp,xh,xl,rts;

	(*funcd)(x1,&fl,&df);
	(*funcd)(x2,&fh,&df);
	if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0))
		nrerror("Root must be bracketed in rtsafe");
	if (fl == 0.0) return x1;
	if (fh == 0.0) return x2;
	if (fl < 0.0) {
		xl=x1;
		xh=x2;
	} else {
		xh=x1;
		xl=x2;
	}
	rts=0.5*(x1+x2);
	dxold=fabs(x2-x1);
	dx=dxold;
	(*funcd)(rts,&f,&df);
	for (j=1;j<=MAXIT2;j++) {
		if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0)
			|| (fabs(2.0*f) > fabs(dxold*df))) {
			dxold=dx;
			dx=0.5*(xh-xl);
			rts=xl+dx;
			if (xl == rts) {printf("iteration %d\n", j); return rts;}
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) {printf("iteration %d\n", j); return rts;}
		}
		if (fabs(dx) < xacc) {printf("iteration %d\n", j); return rts;}

		(*funcd)(rts,&f,&df);
		if (f < 0.0)
			xl=rts;
		else
			xh=rts;
	}
	nrerror("Maximum number of iterations exceeded in rtsafe");
	return 0.0;
}

// driver func
static float fx(float x)
{
        return bessj0(x);
}

// for newton methods
static void funcd(float x,float *fn, float *df)
{
	float ans, d1, d2, e;

	e = sqrt(2000 - (x*x/100)) / 20;
	d1 = (x*exp(-x/200) *sin(e)) /2000 * sqrt(2000 - (x*x/100));
	d2 = (exp(-x/200) *cos(e)) / 200;
	ans = d1 - d2;

	*fn = bessj0(x);
	*df = ans;
}

int main(void)
{
	int i,nb=NBMAX;
	float xacc,root,*xb1,*xb2;

	xb1=vector(1,NBMAX);
	xb2=vector(1,NBMAX);
	zbrak(fx,X1,X2,N,xb1,xb2,&nb);
	
	printf("f(R) = e^0.005R cos(floor(sqrt(2000 - 0.01R^2) * 0.05)) - 0.01\n\n");

	printf("Roots of function f(R) = 0 in the interval [0, 400] when r.e. = 10^-4\n\n");

        printf("Using bisection method:\n");
        for (i=1;i<=nb;i++) {
                xacc=(1.0e-4)*(xb1[i]+xb2[i])/2.0;
                root=rtbis(fx,xb1[i],xb2[i],xacc);
                printf("root: %.6f \nf(R) = %.9f\n",root,fx(root));
        }

        printf("\nUsing linear intersection method:\n");
        for (i=1;i<=nb;i++) {
                xacc=(1.0e-4)*(xb1[i]+xb2[i])/2.0;
                root=rtflsp(fx,xb1[i],xb2[i],xacc);
                printf("root: %.6f \nf(R) = %.9f\n",root,fx(root));
        }
	
	printf("\nUsing Newton with bracketing method:\n");
        for (i=1;i<=nb;i++) {
                xacc=(1.0e-4)*(xb1[i]+xb2[i])/2.0;
                root=rtsafe(funcd,xb1[i],xb2[i],xacc);
                printf("root: %.6f\nf(R) = %.9f\n",root,fx(root));
        }

        printf("\nUsing Muller's method:\n");
        for (i=1;i<=nb;i++) {
                xacc=(1.0e-4)*(xb1[i]+xb2[i])/2.0;
                root=muller(fx,xb1[i],xb2[i],xacc);
                printf("root: %.6f\nf(R) = %.9f\n",root,fx(root));
        }


	printf("\n\nRoots of function f(R) = 0 in the interval [0, 400] when r.e. = 10^-6\n\n");	

	printf("Using bisection method:\n");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtbis(fx,xb1[i],xb2[i],xacc);
		printf("root: %.6f \nf(R) = %.9f\n",root,fx(root));
	}

	printf("\nUsing linear intersection method:\n");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtflsp(fx,xb1[i],xb2[i],xacc);
		printf("root: %.6f \nf(R) = %.9f\n",root,fx(root));
	}
/*
	printf("\nUsing secant method:\n");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsec(fx,xb1[i],xb2[i],xacc);
		printf("root: %.6f \nf(R) = %.9f\n",root,fx(root));
	}

	printf("\nUsing Newton-Raphson method:\n");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtnewt(funcd,xb1[i],xb2[i],xacc);
		printf("root: %.6f\nf(R) = %.9f\n",root,fx(root));
	}
*/
	printf("\nUsing Newton with bracketing method:\n");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(funcd,xb1[i],xb2[i],xacc);
		printf("root: %.6f\nf(R) = %.9f\n",root,fx(root));
	}

	 printf("\nUsing Muller's method:\n");
        for (i=1;i<=nb;i++) {
                xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
                root=muller(fx,xb1[i],xb2[i],xacc);
                printf("root: %.6f\nf(R) = %.9f\n",root,fx(root));
        }

	printf("\n***Root not found using secant, Newton & Muller's method.\n");	
	
	free_vector(xb2,1,NBMAX);
	free_vector(xb1,1,NBMAX);
	return 0;
}
#undef JMAX
#undef MAXIT
#undef MAXIT2
#undef NMAX
#undef MAX
#undef NRANSI
