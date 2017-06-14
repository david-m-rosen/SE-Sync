/*
This is the test file to run the problem defined in DriverElasticCurvesRO.h and DriverElasticCurvesRO.cpp.

---- WH
*/

#ifndef TESTELASTICCURVESRO_H
#define TESTELASTICCURVESRO_H

#include <fstream>
#include "Problems/ElasticCurvesRO/DriverElasticCurvesRO.h"
#include "Others/def.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTELASTICCURVESRO)
int main(void);
#endif

#endif // end of TESTELASTICCURVESRO_H
