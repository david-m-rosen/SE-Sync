/*
This is the test file to run the problem defined in LRMatrixCompletion.h and LRMatrixCompletion.cpp.

---- WH
*/

#ifndef TESTLRMATRIXCOMPLETION_H
#define TESTLRMATRIXCOMPLETION_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "test/DriverMexProb.h"

#include "Problems/LRMatrixCompletion/LRMatrixCompletion.h"
#include "Manifolds/LowRank/LowRank.h"
#include "Manifolds/LowRank/LowRankVariable.h"

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

#if !defined(MATLAB_MEX_FILE) && defined(TESTLRMATRIXCOMPLETION)
int main(void);
#endif

void testLRMatrixCompletion(void);

#endif
