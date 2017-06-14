
#include "test/TestEucQuadratic.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTEUCQUADRATIC)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	genrandseed(0);

	// size of the domain
	integer dim = 10;

	// Generate the matrices in the Euclidean Quadratic problem.
	// Use blas to obtain a positive definite matrix by M = Temp * Temp^T
	double *M = new double[dim * dim];
	double *Temp = new double[dim * dim];
	for (integer i = 0; i < dim * dim; i++)
		Temp[i] = genrandnormal();
	char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
	double one = 1, zero = 0;
	integer N = dim;
	printf("start\n");
	// M = temp * temp^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	dgemm_(transn, transt, &N, &N, &N, &one, Temp, &N, Temp, &N, &zero, M, &N);
	printf("end\n");

	delete[] Temp;

	CheckMemoryDeleted = new std::map<integer *, integer>;

	testEucQuadratic(M, dim);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] M;

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}
#endif

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 2)
    {
        mexErrMsgTxt("The number of arguments should be at least two.\n");
    }
    
	double *M, *X, *Xopt;
	M = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer dim;
	dim = mxGetM(prhs[0]);
    if(mxGetN(prhs[0]) != dim)
    {
        mexErrMsgTxt("The size of matrix is not correct.\n");
    }
    if(mxGetM(prhs[1]) != dim || mxGetN(prhs[1]) != 1)
    {
        mexErrMsgTxt("The size of the initial X is not correct!\n");
    }
    
	printf("dim:%d\n", dim);

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(dim, 1, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testEucQuadratic(M, dim, X, Xopt);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	return;
}

#endif

void testEucQuadratic(double *M, integer dim, double *X, double *Xopt)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	printf("tt:%ud\n", tt);
	tt = 0;
	genrandseed(tt);

	// Obtain an initial iterate
	EucVariable EucX(dim, 1);
	if (X == nullptr)
	{
		EucX.RandInManifold();
	}
	else
	{
		double *EucXptr = EucX.ObtainWriteEntireData();
		for (integer i = 0; i < dim; i++)
			EucXptr[i] = X[i];
	}

	// Define the manifold
	Euclidean Domain(dim);

	// Define the problem
	EucQuadratic Prob(M, dim);
	Prob.SetDomain(&Domain);

	// test RSD
	printf("********************************Check all line search algorithm in RSD*****************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RSD *RSDsolver = new RSD(&Prob, &EucX);
		RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RSDsolver->Debug = FINALRESULT;
		RSDsolver->CheckParams();
		RSDsolver->Run();
		delete RSDsolver;
	}
	// test RNewton
	printf("********************************Check all line search algorithm in RNewton*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RNewton *RNewtonsolver = new RNewton(&Prob, &EucX);
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
		RCG *RCGsolver = new RCG(&Prob, &EucX);
		RCGsolver->RCGmethod = static_cast<RCGmethods> (i);
		RCGsolver->LineSearch_LS = STRONGWOLFE;
		RCGsolver->LS_beta = 0.1;
		RCGsolver->Debug = FINALRESULT;
		RCGsolver->CheckParams();
		RCGsolver->Run();
		delete RCGsolver;
	}

	// test RBroydenFamily
	printf("********************************Check all line search algorithm in RBroydenFamily*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RBroydenFamily *RBroydenFamilysolver = new RBroydenFamily(&Prob, &EucX);
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
		RWRBFGS *RWRBFGSsolver = new RWRBFGS(&Prob, &EucX);
		RWRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RWRBFGSsolver->Debug = FINALRESULT;
		RWRBFGSsolver->CheckParams();
		RWRBFGSsolver->Run();
		delete RWRBFGSsolver;
	}

	// test RBFGS
	printf("********************************Check all line search algorithm in RBFGS*************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RBFGS *RBFGSsolver = new RBFGS(&Prob, &EucX);
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
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &EucX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->Debug = FINALRESULT;
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		delete LRBFGSsolver;
	}

	printf("********************************Check RTRSD*************************************\n");
	RTRSD RTRSDsolver(&Prob, &EucX);
	printf("\n");
	RTRSDsolver.Debug = FINALRESULT;
	RTRSDsolver.CheckParams();
	RTRSDsolver.Run();

	printf("********************************Check RTRNewton*************************************\n");
	RTRNewton RTRNewtonsolver(&Prob, &EucX);
	printf("\n");
	RTRNewtonsolver.Debug = FINALRESULT;
	RTRNewtonsolver.CheckParams();
	RTRNewtonsolver.Run();

	printf("********************************Check RTRSR1*************************************\n");
	RTRSR1 RTRSR1solver(&Prob, &EucX);
	printf("\n");
	RTRSR1solver.Debug = FINALRESULT;
	RTRSR1solver.CheckParams();
	RTRSR1solver.Run();

	printf("********************************Check LRTRSR1*************************************\n");
	LRTRSR1 LRTRSR1solver(&Prob, &EucX);
	printf("\n");
	LRTRSR1solver.Debug = FINALRESULT;
	LRTRSR1solver.CheckParams();
	LRTRSR1solver.Run();

	// Check gradient and Hessian
	Prob.CheckGradHessian(&EucX);
	const Variable *xopt = RTRNewtonsolver.GetXopt();
	Prob.CheckGradHessian(xopt);
    
	if (Xopt != nullptr)
	{
		const double *xoptptr = xopt->ObtainReadData();
		for (integer i = 0; i < dim; i++)
			Xopt[i] = xoptptr[i];
	}
};
