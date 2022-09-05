#include <math.h>
#define MAXIT 30

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
#undef MAXIT
