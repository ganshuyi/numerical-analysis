#include <math.h> // bessj0

#include <stdio.h>
#include <stdlib.h>
#define NRANSI
#include "nr.h"
#include "nrutil.h"

#define N 100
#define NBMAX 20
//intervals
#define X1 1.0
#define X2 10.0

#define JMAX 40 //rtbis

#define MAXIT 30 //rtflsp, rtsec

#define NMAX 20 //rtnewt

#define MAXIT2 100 //rtsafe

#define NR_END 1
#define FREE_ARG char*

// bessel function
float bessj0(float x)
{
	float ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
			+y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
		ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
			+y*(59272.64853+y*(267.8532712+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-0.785398164;
		ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
			+y*(-0.2073370639e-5+y*0.2093887211e-6)));
		ans2 = -0.1562499995e-1+y*(0.1430488765e-3
			+y*(-0.6911147651e-5+y*(0.7621095161e-6
			-y*0.934945152e-7)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
	}
	return ans;
}

float bessj1(float x)
{
	float ax,z;
	double xx,y,ans,ans1,ans2;

	if ((ax=fabs(x)) < 8.0) {
		y=x*x;
		ans1=x*(72362614232.0+y*(-7895059235.0+y*(242396853.1
			+y*(-2972611.439+y*(15704.48260+y*(-30.16036606))))));
		ans2=144725228442.0+y*(2300535178.0+y*(18583304.74
			+y*(99447.43394+y*(376.9991397+y*1.0))));
		ans=ans1/ans2;
	} else {
		z=8.0/ax;
		y=z*z;
		xx=ax-2.356194491;
		ans1=1.0+y*(0.183105e-2+y*(-0.3516396496e-4
			+y*(0.2457520174e-5+y*(-0.240337019e-6))));
		ans2=0.04687499995+y*(-0.2002690873e-3
			+y*(0.8449199096e-5+y*(-0.88228987e-6
			+y*0.105787412e-6)));
		ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
		if (x < 0.0) ans = -ans;
	}
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
		if (fabs(dx) < xacc || fmid == 0.0) return rtb;
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
		if (fabs(del) < xacc || f == 0.0) return rtf;
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
	for (j=1;j<=MAXIT;j++) {
		dx=(xl-rts)*f/(f-fl);
		xl=rts;
		fl=f;
		rts += dx;
		f=(*func)(rts);
		if (fabs(dx) < xacc || f == 0.0) return rts;
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
		if (fabs(dx) < xacc || j > MAXIT) return result;
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
		if (fabs(dx) < xacc) return rtn;
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
			if (xl == rts) return rts;
		} else {
			dxold=dx;
			dx=f/df;
			temp=rts;
			rts -= dx;
			if (temp == rts) return rts;
		}
		if (fabs(dx) < xacc) return rts;
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
	*fn=bessj0(x);
	*df = -bessj1(x);
}

int main(void)
{
	int i,nb=NBMAX;
	float xacc,root,*xb1,*xb2;

	xb1=vector(1,NBMAX);
	xb2=vector(1,NBMAX);
	zbrak(fx,X1,X2,N,xb1,xb2,&nb);
	
	printf("Roots of Bessel function in the interval [1.0, 10.0]\n\n");	

	printf("Using bisection method:\n");
	printf("%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtbis(fx,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %18.9f\n",i,root,fx(root));
	}

	printf("\nUsing linear intersection method:\n");
        printf("%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtflsp(fx,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %18.9f\n",i,root,fx(root));
	}

	printf("\nUsing secant method:\n");
        printf("%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsec(fx,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %18.9f\n",i,root,fx(root));
	}

	printf("\nUsing Newton-Raphson method:\n");
        printf("%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtnewt(funcd,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %18.9f\n",i,root,fx(root));
	}

	printf("\nUsing Newton with bracketing method:\n");
        printf("%21s %15s\n","x","f(x)");
	for (i=1;i<=nb;i++) {
		xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
		root=rtsafe(funcd,xb1[i],xb2[i],xacc);
		printf("root %3d %14.6f %18.9f\n",i,root,fx(root));
	}

	 printf("\nUsing Muller's method:\n");
        printf("%21s %15s\n","x","f(x)");
        for (i=1;i<=nb;i++) {
                xacc=(1.0e-6)*(xb1[i]+xb2[i])/2.0;
                root=muller(fx,xb1[i],xb2[i],xacc);
                printf("root %3d %14.6f %18.9f\n",i,root,fx(root));
        }

	
	
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
