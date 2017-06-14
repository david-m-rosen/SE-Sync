
#include "test/TestSPDTensorDL.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSTIEBROCKETT is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTSPDTENSORDL)

int main(void)
{
	testSPDTensorDL();

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

#endif

void testSPDTensorDL(void)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	printf("tt:%ud\n", tt);//---
	genrandseed(tt);

	/*Randomly generate a point on the SPD manifold*/
	integer dim = 3, num = 2, N = 10;
	double lambdaX = 1;
	SPDTVariable SPDTX(dim, num);
	SPDTX.RandInManifold();
	// Define the manifold
	SPDTensor Domain(dim, num);
	//Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	double *Ls = new double[dim * dim * N + dim * dim + num * N];
	double *tmp = Ls + dim * dim * N;
	double *alpha = tmp + dim * dim;
	integer info;
	for (integer i = 0; i < N; i++)
	{
		for (integer j = 0; j < dim * dim; j++)
			tmp[j] = genrandnormal();

		dgemm_(GLOBAL::N, GLOBAL::T, &dim, &dim, &dim, &GLOBAL::DONE, tmp, &dim, tmp, &dim, &GLOBAL::DZERO, Ls + i * dim * dim, &dim);

		dpotrf_(GLOBAL::L, &dim, Ls + i * dim * dim, &dim, &info);
		if (info != 0)
		{
			printf("Warning: testSPDTensorDL Cholesky decomposition fails with info:%d!\n", info);
		}
		for (integer j = 0; j < dim; j++)
			for (integer k = j + 1; k < dim; k++)
				Ls[j + k * dim + i * dim * dim] = 0;
	}

	for (integer i = 0; i < num * N; i++)
	{
		alpha[i] = genrandreal();
	}

	// Define the problem
	SPDTensorDL Prob(Ls, alpha, dim, N, num, lambdaX);
	/*The domain of the problem is a SPD manifold*/
	Prob.SetDomain(&Domain);

	//Prob.CheckGradHessian(&SPDX);

	/*Output the parameters of the domain manifold*/
	Domain.CheckParams();

	//Prob.CheckGradHessian(&SPDTX);
	//return;

	/*Check the correctness of the manifold operations*/
	//Domain.CheckIntrExtr(&SPDTX);
	//Domain.CheckRetraction(&SPDTX);
	//Domain.CheckDiffRetraction(&SPDTX);
	//Domain.CheckLockingCondition(&SPDTX);
	//Domain.CheckcoTangentVector(&SPDTX);
	//Domain.CheckIsometryofVectorTransport(&SPDTX);
	//Domain.CheckIsometryofInvVectorTransport(&SPDTX);
	//Domain.CheckVecTranComposeInverseVecTran(&SPDTX);
	//Domain.CheckTranHInvTran(&SPDTX);
	//Domain.CheckHaddScaledRank1OPE(&SPDTX);

	// test LRBFGS
	printf("********************************Test Dictionary Learning for SPD Tensor in LRBFGS*************************************\n");
	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &SPDTX);
	LRBFGSsolver->LineSearch_LS = ARMIJO;
	//LRBFGSsolver->Num_pre_funs = 3;
	LRBFGSsolver->Debug = ITERRESULT; //ITERRESULT;// 
	LRBFGSsolver->InitSteptype = ONESTEP;
	LRBFGSsolver->Max_Iteration = 3;
	LRBFGSsolver->Tolerance = 1e-6;
	LRBFGSsolver->OutputGap = 1;
	LRBFGSsolver->LengthSY = 4;
	//LRBFGSsolver->Accuracy = 1e-4;
	LRBFGSsolver->Finalstepsize = 1;
	LRBFGSsolver->CheckParams();
	LRBFGSsolver->Run();
	delete LRBFGSsolver;

	//// test RCG
	//printf("********************************Test Dictionary Learning for SPD Tensor in RCG*************************************\n");
	//RCG *RCGsolver = new RCG(&Prob, &SPDTX);
	//RCGsolver->LineSearch_LS = ARMIJO;
	////RCGsolver->Num_pre_funs = 3;
	//RCGsolver->Debug = FINALRESULT; //ITERRESULT;// 
	//RCGsolver->InitSteptype = BBSTEP;
	//RCGsolver->Max_Iteration = 2000;
	//RCGsolver->Tolerance = 1e-6;
	//RCGsolver->OutputGap = 10;
	//RCGsolver->RCGmethod = POLAK_RIBIERE_MOD;
	//RCGsolver->CheckParams();
	//RCGsolver->Run();
	//delete RCGsolver;

	delete[] Ls;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 6)
	{
		mexErrMsgTxt("The number of arguments should be at least six.\n");
	}
	double *Ls, *alpha, *Xinitial;
	double lambdaX;
	integer dim, N, num, HasHHR;
	Ls = mxGetPr(prhs[0]); /*The data: dim by dim by N*/
	const mwSize *ptrdims = mxGetDimensions(prhs[0]);
	dim = ptrdims[0];
	if (ptrdims[0] != ptrdims[1])
		mexErrMsgTxt("The size of tensor Ls is not correct.\n");
	N = ptrdims[2];
	alpha = mxGetPr(prhs[1]); /*the sparse coding, num by N*/
	num = mxGetM(prhs[1]);
	if (N != mxGetN(prhs[1]))
		mexErrMsgTxt("The size of tensor alpha does not match the size of Ls.\n");
	Xinitial = mxGetPr(prhs[2]); /*initial iterate, dim by dim by num*/
	ptrdims = mxGetDimensions(prhs[2]);
	if (ptrdims[0] != dim || ptrdims[1] != dim)
		mexErrMsgTxt("The size of initial iterate does not match the size of Ls.\n");

	if (mxGetNumberOfDimensions(prhs[2]) > 2 && ptrdims[2] != num)
		mexErrMsgTxt("The size of initial iterate does not match the size of alpha.\n");

	lambdaX = mxGetScalar(prhs[3]);
	HasHHR = static_cast<integer> (mxGetScalar(prhs[4]));

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate from the input
	SPDTVariable *SPDTX = new SPDTVariable(dim, num);
	double *initialX = SPDTX->ObtainWriteEntireData();
	for (integer i = 0; i < dim * dim * num; i++)
	{
		initialX[i] = Xinitial[i];
	}
	SPDTVariable *SPDTsoln = nullptr;
	if (nrhs >= 7)
	{
		double *soln = mxGetPr(prhs[6]);
		SPDTsoln = new SPDTVariable(dim, num);
		double *SPDTsolnptr = SPDTsoln->ObtainWriteEntireData();
		for (integer i = 0; i < num * N; i++)
		{
			SPDTsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	SPDTensor Domain(dim, num);
	// Define the SPDMean problem
	SPDTensorDL Prob(Ls, alpha, dim, N, num, lambdaX);
	Prob.SetDomain(&Domain);
	Domain.SetHasHHR((HasHHR != 0));
	ParseSolverParamsAndOptimizing(prhs[5], &Prob, SPDTX, SPDTsoln, plhs);
	//Domain.CheckParams();
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete SPDTX;
	delete SPDTsoln;
	return;
}

#endif
