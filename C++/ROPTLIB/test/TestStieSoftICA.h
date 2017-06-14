/*
This is the test file to run the problem defined in StieSoftICA.h and StieSoftICA.cpp.

---- WH
*/

#ifndef TESTSTIESOFTICA_H
#define TESTSTIESOFTICA_H

#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>
#include "test/DriverMexProb.h"

#include "Problems/StieSoftICA/StieSoftICA.h"
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

#if !defined(MATLAB_MEX_FILE) && defined(TESTSTIESOFTICA)
int main(void);
#endif

using namespace ROPTLIB;

#endif // end of TESTSTIESOFTICA_H
