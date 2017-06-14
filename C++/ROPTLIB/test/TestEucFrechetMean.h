/*
This is the test file to run the problem defined in EucFrechetMean.h and EucFrechetMean.cpp.

---- WH
*/

#ifndef TESTEUCFRECHETMEAN_H
#define TESTEUCFRECHETMEAN_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Problems/EucFrechetMean/EucFrechetMean.h"

#include "Problems/StieBrockett/StieBrockett.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/Stiefel.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Others/def.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTEUCFRECHETMEAN)
int main(void);
#endif

void testEucFrechetMean(double *Data, double *Weight, integer num, integer dim, double *X = nullptr, double *Xopt = nullptr);

#endif // end of TESTEUCFRECHETMEAN_H
