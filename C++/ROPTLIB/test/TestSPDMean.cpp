
#include "test/TestSPDMean.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSTIEBROCKETT is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDMEAN)

int main(void)
{
	testSPDMean();

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

#endif

void testSPDMean(void)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	genrandseed(tt);

	/*Randomly generate a point on the SPD manifold*/
	integer n = 100, num = 4;
	SPDVariable SPDX(n);
	double *initialX = SPDX.ObtainWriteEntireData();
	for (integer i = 0; i < n; i++)
	{
		for (integer j = 0; j < n; j++)
		{
			initialX[i + j * n] = 0;
		}
		initialX[i + i * n] = 1;
	}

	// Define the manifold
	SPDManifold Domain(n);
	Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	double *Ls = new double[n * n * num + n * n];
	double *tmp = Ls + n * n * num;
	integer info;
	for (integer i = 0; i < num; i++)
	{
		for (integer j = 0; j < n * n; j++)
			tmp[j] = genrandnormal();

		dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, tmp, &n, tmp, &n, &GLOBAL::DZERO, Ls + i * n * n, &n);


		dpotrf_(GLOBAL::L, &n, Ls + i * n * n, &n, &info);
		if (info != 0)
		{
			printf("Warning: TestSPDMean Cholesky decomposition fails with info:%d!\n", info);
		}
		for (integer j = 0; j < n; j++)
			for (integer k = j + 1; k < n; k++)
				Ls[j + k * n + i * n * n] = 0;
	}

	// Define the problem
	SPDMean Prob(Ls, n, num);
	/*The domain of the problem is a SPD manifold*/
	Prob.SetDomain(&Domain);

	//Prob.CheckGradHessian(&SPDX);

	/*Output the parameters of the domain manifold*/
	Domain.CheckParams();

	/*Check the correctness of the manifold operations*/
	//Domain.CheckIntrExtr(&SPDX);
	//Domain.CheckRetraction(&SPDX);
	//Domain.CheckDiffRetraction(&SPDX);
	//Domain.CheckLockingCondition(&SPDX);
	//Domain.CheckcoTangentVector(&SPDX);
	//Domain.CheckIsometryofVectorTransport(&SPDX);
	//Domain.CheckIsometryofInvVectorTransport(&SPDX);
	//Domain.CheckVecTranComposeInverseVecTran(&SPDX);
	//Domain.CheckTranHInvTran(&SPDX);
	//Domain.CheckHaddScaledRank1OPE(&SPDX);

	// test LRBFGS
	printf("********************************Test Geometric mean in LRBFGS*************************************\n");
	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &SPDX);
	LRBFGSsolver->LineSearch_LS = ARMIJO;
	LRBFGSsolver->Debug = ITERRESULT; //ITERRESULT;// 
	LRBFGSsolver->Max_Iteration = 2000;
	LRBFGSsolver->Tolerance = 1e-10;
	LRBFGSsolver->Accuracy = 1e-4;
	LRBFGSsolver->Finalstepsize = 1;
	LRBFGSsolver->CheckParams();
	LRBFGSsolver->Run();
	//// Check gradient and Hessian
	//Prob.CheckGradHessian(&SPDX);
	//const Variable *xopt = LRBFGSsolver->GetXopt();
	//Prob.CheckGradHessian(xopt);

	LRBFGS *LRBFGSsolver2 = new LRBFGS(&Prob, &SPDX, LRBFGSsolver->GetXopt());
	LRBFGSsolver2->LineSearch_LS = ARMIJO;
	LRBFGSsolver2->Debug = ITERRESULT; //ITERRESULT;// 
	LRBFGSsolver2->Max_Iteration = 500;
	LRBFGSsolver2->Tolerance = 1e-5;
	LRBFGSsolver2->Accuracy = 1e-4;
	LRBFGSsolver2->Finalstepsize = 1;
	LRBFGSsolver2->CheckParams();
	LRBFGSsolver2->Run();
	double *dists = LRBFGSsolver2->GetdistSeries();
	ForDebug::Print("Dist:", dists, LRBFGSsolver2->GetlengthSeries());
	delete LRBFGSsolver2;
	delete LRBFGSsolver;

	delete[] Ls;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 4)
	{
		mexErrMsgTxt("The number of arguments should be at least four.\n");
	}
	double *Ls, *X, *Xopt;
	integer n, N, HasHHR;
	Ls = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	n = mxGetM(prhs[1]);
	const mwSize *ptrdims = mxGetDimensions(prhs[0]);
	if (mxGetNumberOfDimensions(prhs[0]) == 2)
		N = 1;
	else
		N = ptrdims[2];
	HasHHR = static_cast<integer> (mxGetScalar(prhs[2]));

	if (ptrdims[1] != n || ptrdims[0] != n)
	{
		mexErrMsgTxt("The size of matrix C is not correct.\n");
	}
	if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != n)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	//testStieSoftICA(Cs, n, p, N, X, Xopt);

	// Obtain an initial iterate from the input
	SPDVariable *SPDX = new SPDVariable(n);
	double *initialX = SPDX->ObtainWriteEntireData();
	for (integer i = 0; i < n * n; i++)
	{
		initialX[i] = X[i];
	}

	SPDVariable *SPDsoln = nullptr;
	if (nrhs >= 5)
	{
		double *soln = mxGetPr(prhs[4]); /*soln: n by p*/
		SPDsoln = new SPDVariable(n);
		double *SPDsolnptr = SPDsoln->ObtainWriteEntireData();
		for (integer i = 0; i < n * n; i++)
		{
			SPDsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	SPDManifold Domain(n);
	// Define the SPDMean problem
	SPDMean Prob(Ls, n, N);
	Prob.SetDomain(&Domain);
	Domain.SetHasHHR((HasHHR != 0));
	ParseSolverParamsAndOptimizing(prhs[3], &Prob, SPDX, SPDsoln, plhs);
	//Domain.CheckParams();
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete SPDX;
	delete SPDsoln;
	return;
}

#endif
