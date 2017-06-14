/*
This is the test file for the bounding box problem defined in OrthBoundingBox.h and OrthBoundingBox.cpp.

---- WH
*/

#ifndef TESTORTHBOUNDINGBOX_H
#define TESTORTHBOUNDINGBOX_H

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
#include "Problems/OrthBoundingBox/OrthBoundingBox.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/OrthGroup/OrthGroupVector.h"
#include "Manifolds/OrthGroup/OrthGroupVariable.h"
#include "Manifolds/OrthGroup/OrthGroup.h"

/*Linesearch based solvers*/
#include "Solvers/SolversLS.h"
#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"
#include "Solvers/RBFGSLPSub.h"
#include "Solvers/LRBFGSLPSub.h"
#include "Solvers/RGS.h"

/*Trust-region based solvers*/
#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

/*The global head file*/
#include "Others/def.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTORTHBOUNDINGBOX)
int main(void);
#endif

/*The main test function*/
void testOrthBoundingBox(double *E, integer d, integer n, double *X = nullptr, double *Xopt = nullptr);

#endif // end of TESTORTHBOUNDINGBOX_H
