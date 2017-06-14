
#include "test/TestGrassRQ.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTGRASSRQ is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTGRASSRQ)

/*Help to check the memory leakage problem. No necesary any more.*/
std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	/*Set the random seed*/
	genrandseed(0);

	// size of the Grassmann manifold
	integer n = 5, p = 3;

	// Generate the matrices in the Rayleigh Quotient problem.
	double *B = new double[n * n];
	/*B is an n by n matrix*/
	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B[i + j * n] = genrandnormal();
			B[j + i * n] = B[i + j * n];
		}
	}

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testGrassRQ(B, n, p);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] B;
#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

/*We don't have to a line search algorithm defined in the solvers. The line seach algorithm can be defined 
here:*/
double LinesearchInput(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *prob, const Solvers *solver)
{ /*For example, simply use one to be the stepsize*/
	return 1;
}

void testGrassRQ(double *B, integer n, integer p, double *X, double *Xopt)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	genrandseed(tt);
	GrassVariable GrassX(n, p);

	if (X == nullptr)
	{/*If X is not defined before, then obtain an initial iterate by taking the Q factor of qr decomposition*/
		GrassX.RandInManifold();
	}
	else
	{/*Otherwise, using the input orthonormal matrix as the initial iterate*/
		double *GrassXptr = GrassX.ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
			GrassXptr[i] = X[i];
	}

	// Define the manifold
	Grassmann Domain(n, p);
	//Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	// Define the Rayleigh Quotient problem
	GrassRQ Prob(B, n, p);
	/*The domain of the problem is a Grassmann manifold*/
	Prob.SetDomain(&Domain);

	/*Output the parameters of the domain manifold*/
	Domain.CheckParams();

	/*Check the correctness of the manifold operations*/
	//Domain.CheckIntrExtr(&GrassX);
	//Domain.CheckRetraction(&GrassX);
	//Domain.CheckDiffRetraction(&GrassX);
	//Domain.CheckLockingCondition(&GrassX);
	//Domain.CheckcoTangentVector(&GrassX);
	//Domain.CheckIsometryofVectorTransport(&GrassX);
	//Domain.CheckIsometryofInvVectorTransport(&GrassX);
	//Domain.CheckVecTranComposeInverseVecTran(&GrassX);
	//Domain.CheckTranHInvTran(&GrassX);
	//Domain.CheckHaddScaledRank1OPE(&GrassX);

	// test RSD
	printf("********************************Check all line search algorithm in RSD*****************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RSD *RSDsolver = new RSD(&Prob, &GrassX);
		RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RSDsolver->Debug = FINALRESULT;
		RSDsolver->Max_Iteration = 2000;
		RSDsolver->CheckParams();
		RSDsolver->Run();
		delete RSDsolver;
	}

	// test RNewton
	printf("********************************Check all line search algorithm in RNewton*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RNewton *RNewtonsolver = new RNewton(&Prob, &GrassX);
		RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RNewtonsolver->Debug = ITERRESULT;
		/*Uncomment following two lines to use the linesearch algorithm defined by the function "LinesearchInput".*/
		//RNewtonsolver->LineSearch_LS = INPUTFUN;
		//RNewtonsolver->LinesearchInput = &LinesearchInput;
		RNewtonsolver->Max_Iteration = 10;
		RNewtonsolver->CheckParams();
		RNewtonsolver->Run();
		delete RNewtonsolver;
	}

	// test RCG
	printf("********************************Check all Formulas in RCG*************************************\n");
	for (integer i = 0; i < RCGMETHODSLENGTH; i++)
	{
		RCG *RCGsolver = new RCG(&Prob, &GrassX);
		RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
		RCGsolver->LineSearch_LS = ARMIJO;
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
		RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &GrassX);
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
		RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &GrassX);
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
		RBFGS *RBFGSsolver = new RBFGS(&Prob, &GrassX);
		RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RBFGSsolver->Debug = FINALRESULT;
		RBFGSsolver->CheckParams();
		RBFGSsolver->Run();
		delete RBFGSsolver;
	}

	// test LRBFGS
	printf("********************************Check all line search algorithm in LRBFGS*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &GrassX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;// 
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		delete LRBFGSsolver;
	}

	// test RTRSD
	printf("********************************Check RTRSD*************************************\n");
	RTRSD RTRSDsolver(&Prob, &GrassX);
	RTRSDsolver.Debug = FINALRESULT;
	RTRSDsolver.Max_Iteration = 5000;
	RTRSDsolver.CheckParams();
	RTRSDsolver.Run();

	// test RTRNewton
	printf("********************************Check RTRNewton*************************************\n");
	RTRNewton RTRNewtonsolver(&Prob, &GrassX);
	RTRNewtonsolver.Debug = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	// test RTRSR1
	printf("********************************Check RTRSR1*************************************\n");
	RTRSR1 RTRSR1solver(&Prob, &GrassX);
	RTRSR1solver.Debug = FINALRESULT;
	RTRSR1solver.CheckParams();
	RTRSR1solver.Run();

	// test LRTRSR1
	printf("********************************Check LRTRSR1*************************************\n");
	LRTRSR1 LRTRSR1solver(&Prob, &GrassX);
	LRTRSR1solver.Debug = FINALRESULT;
	LRTRSR1solver.CheckParams();
	LRTRSR1solver.Run();

	// Check the correctness of gradient and Hessian at the initial iterate
	Prob.CheckGradHessian(&GrassX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	// Check the correctness of gradient and Hessian at the final iterate of RTRNewton method
	Prob.CheckGradHessian(xopt);
    
	//Output the optimizer obtained by RTRNewton method
	if (Xopt != nullptr)
	{
		const double *xoptptr = xopt->ObtainReadData();
		for (integer i = 0; i < n * p; i++)
			Xopt[i] = xoptptr[i];
	}
}

#endif

/*If it is compiled in Matlab, then the following "mexFunction" is used as the entrance.*/
#ifdef MATLAB_MEX_FILE

/*Help to check the memory leakage problem. No necesary any more.*/
std::map<integer *, integer> *CheckMemoryDeleted;

/*This function checks the number and formats of input parameters.
nlhs: the number of output in mxArray format
plhs: the output objects in mxArray format
nrhs: the number of input in mxArray format
prhs: the input objects in mxArray format */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if(nrhs < 4)
	{
		mexErrMsgTxt("The number of arguments should be at least four.\n");
	}
	double *B, *X, *Xopt, *soln;
	B = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer p, n, HasHHR;
	n = mxGetM(prhs[0]);
	p = mxGetN(prhs[1]);

	/*Check the correctness of the inputs*/
	if(mxGetN(prhs[0]) != n)
	{
		mexErrMsgTxt("The size of matrix B is not correct.\n");
	}
	if(mxGetM(prhs[1]) != n)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[2]));

	printf("(n, p):%d,%d\n", n, p);

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	GrassVariable GrassX(n, p);
	double *GrassXptr = GrassX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		GrassXptr[i] = X[i];

	GrassVariable *Grasssoln = nullptr;
	if (nrhs >= 5)
	{
		double *soln = mxGetPr(prhs[4]); /*soln: n by p*/
		Grasssoln = new GrassVariable(n, p);
		double *Grasssolnptr = Grasssoln->ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
		{
			Grasssolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	Grassmann Domain(n, p);

	// Define the Rayleigh Quotient problem
	GrassRQ Prob(B, n, p);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[3], &Prob, &GrassX, Grasssoln, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete Grasssoln;
	return;
}

#endif
