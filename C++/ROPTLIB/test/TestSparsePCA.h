/*
This is the test file to run the problem defined in ObliqueTestSparsePCA.h and ObliqueTestSparsePCA.cpp.

---- WH
*/

#ifndef TESTSPARSEPCA_H
#define TESTSPARSEPCA_H

#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "Problems/ObliqueSparsePCA/ObliqueSparsePCA.h"
#include "Manifolds/Oblique/Oblique.h"
#include "Manifolds/Oblique/ObliqueVariable.h"
#include "Manifolds/Oblique/ObliqueVector.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"
#include "Solvers/LRBFGSLPSub.h"
#include "Solvers/RBFGSLPSub.h"

#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "test/DriverMexProb.h"

#include "Others/def.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPARSEPCA)
int main(void);
#endif

integer GetNumberBetweenC1andC2(const Element *x, double c1, double c2);
void testSparsePCA(double *B, double *Dsq, integer p, integer n, integer r, double mu, double *X = nullptr, double *Xopt = nullptr);

#endif // end of TESTTESTSPARSEPCA_H