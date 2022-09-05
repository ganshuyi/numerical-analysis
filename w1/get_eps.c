#include <stdio.h>
#include <float.h>

void get_eps_f(float eps)
{
	float prev_eps;

	while ((1+eps) != 1)
	{
		prev_eps = eps;
		eps /= 2;
	}
	printf("Machine accuracy of float:  %.17g\n", prev_eps);
}

void get_eps_d(double eps)
{
	double prev_eps;

	while ((1+eps) != 1)
        {
                prev_eps = eps;
                eps /= 2;
        }
        printf("Machine accuracy of double: %.17g\n", prev_eps);
}

int main(void)
{
	get_eps_f(0.5);
	get_eps_d(0.5);


	return 0;
}
