/*
This is the test file to run the problem defined in EucQuadratic.h and EucQuadratic.cpp.

---- WH
*/

#ifndef TESTEUCQUADRATIC_H
#define TESTEUCQUADRATIC_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Problems/EucQuadratic/EucQuadratic.h"

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

#if !defined(MATLAB_MEX_FILE) && defined(TESTEUCQUADRATIC)
int main(void);
#endif

void testEucQuadratic(double *M, integer dim, double *X = nullptr, double *Xopt = nullptr);

#endif // end of TESTEUCQUADRATIC_H
