
#ifndef RANDGEN_CPP
#define RANDGEN_CPP

#include "Others/randgen.h"

void genrandseed(unsigned int s)
{
    srand(s);
}

double genrandreal(void)
{
    return static_cast<double> (rand()) / RAND_MAX;
}

double genrandnormal(void)
{
	static double rand1, rand2;
	double tmp = genrandreal();
	while (tmp == 1.0)
		tmp = genrandreal();
    rand1 = -2 * log(1.0 - tmp);
    rand2 = (1.0 - genrandreal()) * 6.2831853071795864769252866;
	return sqrt(rand1) * cos(rand2);
}

#endif
