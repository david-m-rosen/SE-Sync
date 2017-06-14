#include "test/TestLRMatrixCompletion.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTLRMATRIXCOMPLETION)

std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	_CrtSetDbgFlag(_CRTDBG_LEAK_CHECK_DF); /*This can detect the memory leakage for global variables!!*/
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	//_CrtSetBreakAlloc(153);
	long seed = static_cast<long> (time(NULL));
	seed = 0;
	printf("seed:%ld\n", seed);
	genrandseed(seed);
	
	CheckMemoryDeleted = new std::map<integer *, integer>;
    
	testLRMatrixCompletion();
	
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	
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
	if (nrhs < 5)
	{
		mexErrMsgTxt("The number of arguments should be at least five.\n");
	}
	double *A, *X, *Xopt;
	A = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer m, n, HasHHR, nzmax, r;
	size_t *ir, *jc;
	nzmax = mxGetNzmax(prhs[0]);
	ir = mxGetIr(prhs[0]);
	jc = mxGetJc(prhs[0]);
	m = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);
	r = static_cast<integer> (mxGetScalar(prhs[2]));

	/*Check the correctness of the inputs*/
	if (mxGetM(prhs[1]) != (m + n + r) * r || mxGetN(prhs[1]) != 1)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[3]));

	printf("m:%d, n:%d\n", m, n);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	LowRankVariable LRX(m, n, r);
	double *LRXptr = LRX.ObtainWriteEntireData();
	for (integer i = 0; i < (m + n + r) * r; i++)
		LRXptr[i] = X[i];

	LowRankVariable *LRsoln = nullptr;
	if (nrhs >= 6)
	{
		double *soln = mxGetPr(prhs[5]); /*soln: n by p*/
		LRsoln = new LowRankVariable(m, n, r);
		double *LRsolnptr = LRsoln->ObtainWriteEntireData();
		for (integer i = 0; i < (m + n + r) * r; i++)
		{
			LRsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	LowRank Domain(m, n, r);
	
	// Define the matrix completion problem
	integer *inir = new integer[2 * nzmax];
	integer *injc = inir + nzmax;

	for (integer i = 0; i < n; i++)
	{
		for (unsigned long long j = jc[i]; j < jc[i + 1]; j++)
		{
			/*row: ir[j], column: i, entry: A[j]*/
			inir[j] = ir[j];
			injc[j] = i;
		}
	}

	LRMatrixCompletion Prob(inir, injc, A, nzmax, m, n, r);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[4], &Prob, &LRX, LRsoln, plhs);

	delete[] inir;

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete LRsoln;
	return;
}

#endif

void testLRMatrixCompletion(void)
{
	integer m = 5, n = 4, r = 2;

	LowRank Domain(m, n, r);
	Domain.SetHasHHR(true);
	LowRankVariable InitialX(m, n, r);
	InitialX.RandInManifold();
	Domain.CheckParams();
	//Domain.CheckIntrExtr(&InitialX);
	//Domain.CheckRetraction(&InitialX);
	//Domain.CheckcoTangentVector(&InitialX);
	//Domain.CheckDiffRetraction(&InitialX, false);
	//Domain.CheckIsometryofVectorTransport(&InitialX);

	//Domain.CheckLockingCondition(&InitialX);
	//Domain.CheckIsometryofInvVectorTransport(&InitialX);
	//Domain.CheckVecTranComposeInverseVecTran(&InitialX);
	//Domain.CheckTranHInvTran(&InitialX);
	//return;

	//InitialX.Print("initialX:");

	// Generate the matrices in the matrix completion approximation problem.
	integer dim = (m + n - r) * r;
	integer nz = 1 * dim;
	integer *ir = new integer[nz * 2];
	integer *jc = ir + nz;

	integer *tmpforidx = new integer[m * n];
	for (integer i = 0; i < m * n; i++)
		tmpforidx[i] = i;
	/*nz number of indices*/
	integer idx = 0, itmp;
	for (integer i = 0; i < nz; i++)
	{
		/*idx is an integer in [0, m - i - 1]*/
		idx = static_cast<integer> ((m * n - i) * genrandreal());
		while (idx >= m * n - i)
			idx = static_cast<integer> ((m * n - i) * genrandreal());
		/*the chosen idx is put at the end of the array*/
		itmp = tmpforidx[m * n - i - 1];
		tmpforidx[m * n - i - 1] = tmpforidx[idx];
		tmpforidx[idx] = itmp;
	}
	for (integer i = 0; i < nz; i++)
	{
		/*tmpforidx[nz - 1 - i]*/
		ir[i] = static_cast<integer> (tmpforidx[nz - 1 - i] / n);
		jc[i] = tmpforidx[nz - 1 - i] - n * ir[i];
	}
	delete[] tmpforidx;

	integer mn = m * n, mr = m * r, nr = n * r;
	double *A_U = new double[mr];
	double *A_V = new double[nr];
	for (integer i = 0; i < m * r; i++)
	{
		A_U[i] = genrandnormal();
	}
	for (integer i = 0; i < n * r; i++)
	{
		A_V[i] = genrandnormal();
	}
	double *V = new double[nz];
	for (integer i = 0; i < nz; i++)
	{
		V[i] = 0;
		for (integer j = 0; j < r; j++)
		{
			V[i] += A_U[ir[i] + j * m] * A_V[jc[i] + j * n];
		}
	}
	delete[]A_U;
	delete[]A_V;

	LRMatrixCompletion Prob(ir, jc, V, nz, m, n, r);
	Prob.SetDomain(&Domain);

	//Prob.f(&InitialX);
	//Vector *gf = Domain.GetEMPTYINTR()->ConstructEmpty();//---
	//Prob.Grad(&InitialX, gf);//---
	//delete gf;

	//Prob.CheckGradHessian(&InitialX);//--

	//RSD *RSDsolver = new RSD(&Prob, &InitialX);
	RTRNewton *RSDsolver = new RTRNewton(&Prob, &InitialX);
	//LRBFGS *RSDsolver = new LRBFGS(&Prob, &InitialX);
	//->LineSearch_LS = ARMIJO;
	//RSDsolver->LS_beta = 0.01;
	//RSDsolver->RCGmethod = DAI_YUAN;
	RSDsolver->Debug = ITERRESULT;
	RSDsolver->OutputGap = 100;
	RSDsolver->Max_Iteration = 500;
	RSDsolver->CheckParams();
	//RSDsolver->Accuracy = 1e-6;
	RSDsolver->Tolerance = 1e-10;
	RSDsolver->Run();
	//Prob.CheckGradHessian(&InitialX);//--
	//Prob.CheckGradHessian(RSDsolver->GetXopt());//--

	RTRNewton *LRBFGSsolver = new RTRNewton(&Prob, &InitialX, RSDsolver->GetXopt());
	LRBFGSsolver->Debug = ITERRESULT;
	LRBFGSsolver->OutputGap = 10;
	LRBFGSsolver->Tolerance = 1e-8;
	LRBFGSsolver->Run();
	for (integer i = 0; i < LRBFGSsolver->GetlengthSeries(); i++)
	{
		std::cout << i << ":" << LRBFGSsolver->GetdistSeries()[i] << std::endl;
	}

	delete RSDsolver;

	delete[] V;
	delete[] ir;

	
};


