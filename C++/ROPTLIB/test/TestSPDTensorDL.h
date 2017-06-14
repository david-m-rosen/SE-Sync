/*
This is the test file for the Dictionary Learning problem on the manifold defined in SPDTensor.h and SPDTensor.cpp.
See details in Section IV.A in [CS15].
	[CS15] Anoop Cherian and Suvrit Sra. "Riemannian Dictionary Learning and Sparse Coding for Positive
	Definite Matrices"

---- WH
*/

#ifndef TESTSPDTENSORDL_H
#define TESTSPDTENSORDL_H

/*Output to console*/
#include <iostream>

/*Generate random number*/
#include "Others/randgen.h"

/*Computational time*/
#include <ctime>

/*If this test file is called from Matlab, then functions in DriverMexProb.h are used.*/
#include "test/DriverMexProb.h"

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/SPDTensorDL/SPDTensorDL.h"

/*Manifold related classes*/
#include "Manifolds/ProductManifold.h"
#include "Manifolds/SPDTensor/SPDTVector.h"
#include "Manifolds/SPDTensor/SPDTVariable.h"
#include "Manifolds/SPDTensor/SPDTensor.h"

/*Linesearch based solvers*/
#include "Solvers/SolversLS.h"
#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

/*Trust-region based solvers*/
#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

/*The global head file*/
#include "Others/def.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDTENSORDL)
int main(void);
#endif

/*The main test function*/
void testSPDTensorDL(void);

#endif // end of TESTSPDTENSORDL_H
