
#include "test/TestEucPosSpCd.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSTIEBROCKETT is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTEUCPOSSPCD)

int main(void)
{
	testEucPosSpCd();

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

#endif

void testEucPosSpCd(void)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 1457100141;
	printf("tt:%ud\n", tt);//---
	genrandseed(tt);

	/*Randomly generate a point on the EucPositive convex set*/
	integer dim = 3, num = 5, N = 1;
	EucPosVariable EPVinit(num, N), EPVtrue(num, N);
	EPVinit.RandInManifold();
	EPVtrue.RandInManifold();
	double *EPVtrueptr = EPVtrue.ObtainWritePartialData();
	EPVtrueptr[0] = 0;
	EPVtrueptr[3] = 0;
	//EPVtrueptr[5] = 0;
	//EPVtrueptr[10] = 0;
	//EPVtrueptr[17] = 0;
	// Define the domain
	EucPositive Domain(num, N);
	//Domain.SetHasHHR(true); /*set whether the idea in [HGA2015, Section 4.3] is used or not*/

	double *Bs = new double[dim * dim * num + dim * dim + dim * dim * N];
	double *tmp = Bs + dim * dim * num;
	double *Ls = tmp + dim * dim;
	integer info;
	for (integer i = 0; i < num; i++)
	{
		for (integer j = 0; j < dim * dim; j++)
			tmp[j] = genrandnormal();

		dgemm_(GLOBAL::N, GLOBAL::T, &dim, &dim, &dim, &GLOBAL::DONE, tmp, &dim, tmp, &dim, &GLOBAL::DZERO, Bs + i * dim * dim, &dim);
	}

	/*Atom = B * EPVtrue */
	integer dd = dim * dim;
	dgemm_(GLOBAL::N, GLOBAL::N, &dd, &N, &num, &GLOBAL::DONE, Bs, &dd, EPVtrue.ObtainWritePartialData(), &num, 
		&GLOBAL::DZERO, Ls, &dd);
	for (integer i = 0; i < N; i++)
	{
		/*Atom is represented by its Cholesky decomposition matrix:L*/
		dpotrf_(GLOBAL::L, &dim, Ls + i * dim * dim, &dim, &info);
		if (info != 0)
		{
			printf("Warning: testEucPosSpCd Cholesky decomposition fails with info:%d!\n", info);
		}
		for (integer j = 0; j < dim; j++)
			for (integer k = j + 1; k < dim; k++)
				Ls[j + k * dim + i * dim * dim] = 0;
	}

	// Define the problem
	EucPosSpCd Prob(Ls, Bs, 0, dim, num, N);

	/*The domain of the problem is a SPD manifold*/
	Prob.SetDomain(&Domain);

	/*Output the parameters of the domain manifold*/
	Domain.CheckParams();

	double *EPVinitptr = EPVinit.ObtainWritePartialData();

	//Prob.CheckGradHessian(&EPVinit);

	///*Check the correctness of the manifold operations*/
	//Domain.CheckIntrExtr(&EPVinit);
	//Domain.CheckRetraction(&EPVinit);
	//Domain.CheckDiffRetraction(&EPVinit);
	//Domain.CheckLockingCondition(&EPVinit);
	//Domain.CheckcoTangentVector(&EPVinit);
	//Domain.CheckIsometryofVectorTransport(&EPVinit);
	//Domain.CheckIsometryofInvVectorTransport(&EPVinit);
	//Domain.CheckVecTranComposeInverseVecTran(&EPVinit);
	//Domain.CheckTranHInvTran(&EPVinit);
	//Domain.CheckHaddScaledRank1OPE(&EPVinit);

	// test LRBFGS
	printf("********************************Test Dictionary Learning for SPD Tensor in LRBFGS*************************************\n");
	LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &EPVinit); //EPVinit EPVtrue
	LRBFGSsolver->LineSearch_LS = ARMIJO;
	//LRBFGSsolver->Num_pre_funs = 3;
	LRBFGSsolver->Debug = ITERRESULT; //ITERRESULT;// 
	LRBFGSsolver->InitSteptype = ONESTEP;
	LRBFGSsolver->Max_Iteration = 500;
	LRBFGSsolver->Tolerance = 1e-6;
	LRBFGSsolver->OutputGap = 1;
	LRBFGSsolver->LengthSY = 4;
	//LRBFGSsolver->Accuracy = 1e-4;
	LRBFGSsolver->Finalstepsize = 1;
	LRBFGSsolver->CheckParams();
	LRBFGSsolver->Run();
	EPVtrue.Print("EPVtrue:");
	LRBFGSsolver->GetXopt()->Print("LRBFGS solution:");
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

	delete[] Bs;
};

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 6)
	{
		mexErrMsgTxt("The number of arguments should be at least six.\n");
	}
	double *Ls, *Bs, *Xinitial;
	double lambdaa;
	integer dim, num, N, HasHHR;
	Ls = mxGetPr(prhs[0]); /*dim by dim by N*/
	const mwSize *ptrdims = mxGetDimensions(prhs[0]);
	dim = ptrdims[0];
	N = (mxGetNumberOfDimensions(prhs[0]) > 2) ? ptrdims[2] : 1;
	if (dim != ptrdims[1])
		mexErrMsgTxt("The size of L is not correct.\n");

	Bs = mxGetPr(prhs[1]); /*Bs: dim by dim by num*/
	ptrdims = mxGetDimensions(prhs[1]);
	num = ptrdims[2];
	if (dim != ptrdims[0] || dim != ptrdims[1])
		mexErrMsgTxt("The size of Bs is not correct.\n");

	Xinitial = mxGetPr(prhs[2]); /*Xinitial: num by N*/
	if (num != mxGetM(prhs[2]) || N != mxGetN(prhs[2]))
		mexErrMsgTxt("The size of Xinitial is not correct.\n");

	lambdaa = static_cast<double> (mxGetScalar(prhs[3]));
	HasHHR = static_cast<integer> (mxGetScalar(prhs[4]));

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate from the input
	EucPosVariable *EPVinit = new EucPosVariable(num, N);
	double *initialX = EPVinit->ObtainWriteEntireData();
	for (integer i = 0; i < num * N; i++)
	{
		initialX[i] = Xinitial[i];
	}
	EucPosVariable *EPVsoln = nullptr;
	if (nrhs >= 7)
	{
		double *soln = mxGetPr(prhs[6]); /*soln: num by N*/
		EPVsoln = new EucPosVariable(num, N);
		double *EPVsolnptr = EPVsoln->ObtainWriteEntireData();
		for (integer i = 0; i < num * N; i++)
		{
			EPVsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	EucPositive Domain(num, N);
	// Define the SPDMean problem
	EucPosSpCd Prob(Ls, Bs, lambdaa, dim, num, N);
	Prob.SetDomain(&Domain);
	Domain.SetHasHHR((HasHHR != 0));
	ParseSolverParamsAndOptimizing(prhs[5], &Prob, EPVinit, EPVsoln, plhs);
	//Domain.CheckParams();
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete EPVinit;
	delete EPVsoln;
	return;
}

#endif
