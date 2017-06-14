/*
This is the test file for the Reyleigh Quotient problem defined in GrassRQ.h and GrassRQ.cpp.

---- WH
*/

#ifndef TESTGRASSRQ_H
#define TESTGRASSRQ_H

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
#include "Problems/GrassRQ/GrassRQ.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/Grassmann/GrassVector.h"
#include "Manifolds/Grassmann/GrassVariable.h"
#include "Manifolds/Grassmann/Grassmann.h"

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

#if !defined(MATLAB_MEX_FILE) && defined(TESTGRASSRQ)
int main(void);
#endif

/*The main test function*/
void testGrassRQ(double *B, integer n, integer p, double *X = nullptr, double *Xopt = nullptr);

#endif // end of TESTGRASSRQ_H
