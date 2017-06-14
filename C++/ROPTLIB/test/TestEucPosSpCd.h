/*
This is the test file for the Sparse coding problem on the manifold defined in EucPositive.h and EucPositive.cpp.
See details in Section IV.B in [CS15].
	[CS15] Anoop Cherian and Suvrit Sra. "Riemannian Dictionary Learning and Sparse Coding for Positive
	Definite Matrices"

---- WH
*/

#ifndef TESTEUCPOSSPCD_H
#define TESTEUCPOSSPCD_H

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
#include "Problems/EucPosSpCd/EucPosSpCd.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/EucPositive/EucPositive.h"
#include "Manifolds/EucPositive/EucPosVariable.h"
#include "Manifolds/EucPositive/EucPosVector.h"

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

#if !defined(MATLAB_MEX_FILE) && defined(TESTEUCPOSSPCD)
int main(void);
#endif

/*The main test function*/
void testEucPosSpCd(void);

#endif // end of TESTEUCPOSSPCD_H
