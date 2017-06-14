
/*wrapper for the C++ default random generator
 * it can be replaced by any other random generator if necessary
 *
 *
 * --by Wen Huang*/

#ifndef RANDGEN_H
#define RANDGEN_H

//#include <cmath>
#include <random>
#include <stdio.h>
#include <iostream>

/* initializes  a seed */
void genrandseed(unsigned int s);

/* generates a random number on [0,1]-real-interval */
double genrandreal();

/* generate a random number following the standard normal distribution*/
double genrandnormal();

#endif
