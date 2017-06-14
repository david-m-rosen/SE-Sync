/*
This is the test file for the Brocokett problem defined in StieSumBrockett.h and StieSumBrockett.cpp.

---- WH
*/

#ifndef TESTPRODUCTEXAMPLE_CPP
#define TESTPRODUCTEXAMPLE_CPP

/*Output to console*/
#include <iostream>
/*Generate random number*/
#include "Others/randgen.h"
/*Computational time*/
#include <ctime>

/*Problem related classes*/
#include "Problems/Problem.h"
#include "Problems/StieSumBrockett/StieSumBrockett.h"

/*Manifold related classes*/
#include "Manifolds/Manifold.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"

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

/*If the file is compiled in Matlab or TESTPRODUCTEXAMPLE is defined in def.h file, then the following
main() function will be called. */
#if defined(MATLAB_MEX_FILE) || defined(TESTPRODUCTEXAMPLE)

void testProductExample();

int main(void)
{
	testProductExample();

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

void testProductExample()
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	genrandseed(tt);

	// size of the Stiefel manifold
	integer n = 12, p = 8, m = 6, q = 2;

	// Generate the matrices in the Brockett problem.
	double *B1 = new double[n * n * 2 + p * 2 + m * m + q];
	double *B2 = B1 + n * n;
	double *B3 = B2 + n * n;
	double *D1 = B3 + m * m;
	double *D2 = D1 + p;
	double *D3 = D2 + p;

	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B1[i + j * n] = genrandnormal();
			B1[j + i * n] = B1[i + j * n];

			B2[i + j * n] = genrandnormal();
			B2[j + i * n] = B2[i + j * n];
		}
	}
	for (integer i = 0; i < m; i++)
	{
		for (integer j = i; j < m; j++)
		{
			B3[i + j * m] = genrandnormal();
			B3[j + i * m] = B3[i + j * m];
		}
	}
	for (integer i = 0; i < p; i++)
	{
		D1[i] = static_cast<double> (i + 1);
		D2[i] = D1[i];
	}
	for (integer i = 0; i < q; i++)
	{
		D3[i] = static_cast<double> (i + 1);
	}

	// number of manifolds in product of manifold
	integer numofmanis = 2; // two kinds of manifolds
	integer numofmani1 = 2; // the first one has two
	integer numofmani2 = 1; // the second one has one

	// Obtain an initial iterate
	StieVariable StieX1(n, p);
	StieVariable StieX2(m, q);
	ProductElement ProdX(numofmanis, &StieX1, numofmani1, &StieX2, numofmani2);
	ProdX.RandInManifold();

	// Define the Stiefel manifold
	Stiefel mani1(n, p);
	Stiefel mani2(m, q);
	ProductManifold Domain(numofmanis, &mani1, numofmani1, &mani2, numofmani2);

	// Define the Brockett problem
	StieSumBrockett Prob(B1, D1, B2, D2, B3, D3, n, p, m, q);

	// Set the domain of the problem to be the Stiefel manifold
	Prob.SetDomain(&Domain);

	// output the parameters of the manifold of domain
	Domain.CheckParams();

	// test RSD
	printf("********************************Check all line search algorithm in RSD*****************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RSD *RSDsolver = new RSD(&Prob, &ProdX);
		RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RSDsolver->Debug = FINALRESULT;// FINALRESULT;
		RSDsolver->CheckParams();
		RSDsolver->Run();
		delete RSDsolver;
	}

	// test RNewton
	printf("********************************Check all line search algorithm in RNewton*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RNewton *RNewtonsolver = new RNewton(&Prob, &ProdX);
		RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RNewtonsolver->Debug = FINALRESULT;
		RNewtonsolver->CheckParams();
		RNewtonsolver->Run();
		delete RNewtonsolver;
	}

	// test RCG
	printf("********************************Check all Formulas in RCG*************************************\n");
	for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	{
		RCG *RCGsolver = new RCG(&Prob, &ProdX);
		RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
		RCGsolver->LineSearch_LS = STRONGWOLFE;
		RCGsolver->LS_beta = 0.1;
		RCGsolver->Debug = FINALRESULT;
		RCGsolver->CheckParams();
		RCGsolver->Run();
		delete RCGsolver;
	}

	// test RBroydenFamily
	printf("********************************Check all Formulas in RCG*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &ProdX);
		RBroydenFamilysolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RBroydenFamilysolver->Debug = FINALRESULT;
		RBroydenFamilysolver->CheckParams();
		RBroydenFamilysolver->Run();
		delete RBroydenFamilysolver;
	}

	// test RWRBFGS
	printf("********************************Check all line search algorithm in RWRBFGS*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &ProdX);
		RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RWRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
		RWRBFGSsolver->CheckParams();
		RWRBFGSsolver->Run();
		delete RWRBFGSsolver;
	}

	// test RBFGS
	printf("********************************Check all line search algorithm in RBFGS*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RBFGS *RBFGSsolver = new RBFGS(&Prob, &ProdX);
		RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RBFGSsolver->Debug = FINALRESULT;
		RBFGSsolver->CheckParams();
		RBFGSsolver->Run();
		delete RBFGSsolver;
	}

	// test LRBFGS
	printf("********************************Check all line search algorithm in LRBFGS*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)//LSALGOLENGTH
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &ProdX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;//
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		delete LRBFGSsolver;
	}

	// test RTRSD
	printf("********************************Check RTRSD*************************************\n");
	RTRSD RTRSDsolver(&Prob, &ProdX);
	RTRSDsolver.Debug = FINALRESULT;
	RTRSDsolver.CheckParams();
	RTRSDsolver.Run();

	// test RTRNewton
	printf("********************************Check RTRNewton*************************************\n");
	RTRNewton RTRNewtonsolver(&Prob, &ProdX);
	RTRNewtonsolver.Debug = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	// test RTRSR1
	printf("********************************Check RTRSR1*************************************\n");
	RTRSR1 RTRSR1solver(&Prob, &ProdX);
	RTRSR1solver.Debug = FINALRESULT;
	RTRSR1solver.CheckParams();
	RTRSR1solver.Run();

	// test LRTRSR1
	printf("********************************Check LRTRSR1*************************************\n");
	LRTRSR1 LRTRSR1solver(&Prob, &ProdX);
	LRTRSR1solver.Debug = FINALRESULT;
	LRTRSR1solver.CheckParams();
	LRTRSR1solver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&ProdX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);

	delete[] B1;

	return;
}

#endif // end of TESTPRODUCTEXAMPLE
#endif // end of TESTPRODUCTEXAMPLE_CPP
