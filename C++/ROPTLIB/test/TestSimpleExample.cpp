/*
This is the test file for the Brocokett problem defined in StieBrockett.h and StieBrockett.cpp.

---- WH
*/

#ifndef TESTSIMPLEEXAMPLE_CPP
#define TESTSIMPLEEXAMPLE_CPP

/*Output to console*/
#include <iostream>
/*Generate random number*/
#include "Others/randgen.h"
/*Computational time*/
#include <ctime>

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/StieBrockett/StieBrockett.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/Stiefel.h"

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

/*If the file is compiled in Matlab , then the following mexFunction() function is the entrance. */
#ifdef MATLAB_MEX_FILE

int main(void);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    main();
	return;
}
#endif

/*If the file is compiled in Matlab or TESTSIMPLEEXAMPLE is defined in def.h file, then the following
main() function will be called. */
#if defined(MATLAB_MEX_FILE) || defined(TESTSIMPLEEXAMPLE)
int main(void)
{
	// choose a random seed based on current clock
	unsigned tt = (unsigned)time(NULL);
	genrandseed(0);

	// size of the Stiefel manifold
	integer n = 12, p = 8;

	// Generate the matrices in the Brockett problem.
	double *B = new double[n * n + p];
	double *D = B + n * n;
	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B[i + j * n] = genrandnormal();
			B[j + i * n] = B[i + j * n];
		}
	}
	for (integer i = 0; i < p; i++)
		D[i] = static_cast<double> (i + 1);

	// Obtain an initial iterate
	StieVariable StieX(n, p);
	StieX.RandInManifold();

	// Define the Stiefel manifold
	Stiefel Domain(n, p);

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);

	// Set the domain of the problem to be the Stiefel manifold
	Prob.SetDomain(&Domain);

	// output the parameters of the manifold of domain
	Domain.CheckParams();


	//// test RSD
	//printf("********************************Check all line search algorithm in RSD*****************************************\n");
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RSD *RSDsolver = new RSD(&Prob, &StieX);
	//	RSDsolver->Debug = FINALRESULT;
	//	RSDsolver->CheckParams();
	//	RSDsolver->Run();
	//	delete RSDsolver;
	//}

	//// test RNewton
	//printf("********************************Check all line search algorithm in RNewton*************************************\n");
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RNewton *RNewtonsolver = new RNewton(&Prob, &StieX);
	//	RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RNewtonsolver->Debug = FINALRESULT;
	//	RNewtonsolver->CheckParams();
	//	RNewtonsolver->Run();
	//	delete RNewtonsolver;
	//}

	//// test RCG
	//printf("********************************Check all Formulas in RCG*************************************\n");
	//for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	//{
	//	RCG *RCGsolver = new RCG(&Prob, &StieX);
	//	RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
	//	RCGsolver->LineSearch_LS = STRONGWOLFE;
	//	RCGsolver->LS_beta = 0.1;
	//	RCGsolver->Debug = FINALRESULT;
	//	RCGsolver->CheckParams();
	//	RCGsolver->Run();
	//	delete RCGsolver;
	//}

	//// test RBroydenFamily
	//printf("********************************Check all Formulas in RBroydenFamily*************************************\n");
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &StieX);
	//	RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RBroydenFamilysolver->Debug = FINALRESULT;
	//	RBroydenFamilysolver->CheckParams();
	//	RBroydenFamilysolver->Run();
	//	delete RBroydenFamilysolver;
	//}

	//// test RWRBFGS
	//printf("********************************Check all line search algorithm in RWRBFGS*************************************\n");
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &StieX);
	//	RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RWRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
	//	RWRBFGSsolver->CheckParams();
	//	RWRBFGSsolver->Run();
	//	delete RWRBFGSsolver;
	//}

	//// test RBFGS
	//printf("********************************Check all line search algorithm in RBFGS*************************************\n");
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RBFGS *RBFGSsolver = new RBFGS(&Prob, &StieX);
	//	RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RBFGSsolver->Debug = FINALRESULT;
	//	RBFGSsolver->CheckParams();
	//	RBFGSsolver->Run();
	//	delete RBFGSsolver;
	//}

	//// test LRBFGS
	//printf("********************************Check all line search algorithm in LRBFGS*************************************\n");
	//for (integer i = 0; i < INPUTFUN; i++)//LSALGOLENGTH
	//{
	//	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &StieX);
	//	LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
	//	LRBFGSsolver->CheckParams();
	//	LRBFGSsolver->Run();
	//	delete LRBFGSsolver;
	//}

	//// test RTRSD
	//printf("********************************Check RTRSD*************************************\n");
	//RTRSD RTRSDsolver(&Prob, &StieX);
	//RTRSDsolver.Debug = FINALRESULT;
	//RTRSDsolver.CheckParams();
	//RTRSDsolver.Run();

	// test RTRNewton
	printf("********************************Check RTRNewton*************************************\n");
	RTRNewton RTRNewtonsolver(&Prob, &StieX);
	RTRNewtonsolver.Debug = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	//// test RTRSR1
	//printf("********************************Check RTRSR1*************************************\n");
	//RTRSR1 RTRSR1solver(&Prob, &StieX);
	//RTRSR1solver.Debug = FINALRESULT;
	//RTRSR1solver.CheckParams();
	//RTRSR1solver.Run();

	//// test LRTRSR1
	//printf("********************************Check LRTRSR1*************************************\n");
	//LRTRSR1 LRTRSR1solver(&Prob, &StieX);
	//LRTRSR1solver.Debug = FINALRESULT;
	//LRTRSR1solver.CheckParams();
	//LRTRSR1solver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&StieX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);

	delete[] B;

	return 0;
}

#endif // end of TESTSIMPLEEXAMPLE
#endif // end of TESTSIMPLEEXAMPLE_CPP
