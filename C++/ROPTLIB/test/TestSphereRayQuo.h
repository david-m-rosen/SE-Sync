/*
This is the test file to run the problem defined in StieBrockett.h and StieBrockett.cpp with p = 1, 
i.e., \min_{x \in S^{n - 1}) x^T B x.

---- WH
*/

#ifndef TESTSPHERERAYQUO_H
#define TESTSPHERERAYQUO_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

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

#include "Manifolds/Sphere/Sphere.h"
#include "Manifolds/Sphere/SphereVariable.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPHERERAYQUO)
int main(void);
#endif

void testSphereRayQuo(double *B, double *D, integer n, double *X = nullptr, double *Xopt = nullptr);

#endif // end of TESTSPHERERAYQUO_H