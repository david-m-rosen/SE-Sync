/*
This is the test file to check the correctenss of the manifold defined in CpxNSTQOrth.h and CpxNSTQOrth.cpp.

---- WH
*/

#ifndef TESTCSO_H
#define TESTCSO_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>
#include "Test/DriverMexProb.h"

#include "Manifolds/CpxNStQOrth/CSOVector.h"
#include "Manifolds/CpxNStQOrth/CSOVariable.h"
#include "Manifolds/CpxNStQOrth/CpxNStQOrth.h"

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

#if !defined(MATLAB_MEX_FILE) && defined(TESTCSO)
int main(void);
#endif

void testCSO();

#endif // end of TESTCSO_H
