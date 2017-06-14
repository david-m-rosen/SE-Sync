
#include "test/TestStieSoftICA.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTSTIESOFTICA)

std::map<integer *, integer> *CheckMemoryDeleted;

void testStieSoftICA(double *Cs, integer n, integer p, integer N, double *X = nullptr, double *Xopt = nullptr);

int main(void)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 1466602277;
	printf("seed:%ul\n", tt);
	genrandseed(tt);

	// size of the Stiefel manifold
	integer n = 64, p = 32;
	// number of covariance matrices
	integer N = 32;

	// Generate the matrices in the joint diagonalization (JD) problem (Soft ICA problem).
	// Use the same approach as the experiments in paper "A Riemannian symmetric rank-one trust-region method".
	double *Cs = new double[n * n * N];
	for (integer i = 0; i < N; i++)
	{
		for (integer j = 0; j < n; j++)
		{
			for (integer k = 0; k < n; k++)
			{
				Cs[i * n * n + j * n + k] = genrandnormal();
			}
		}
		for (integer j = 0; j < n; j++)
		{
			for (integer k = j; k < n; k++)
			{
				Cs[i * n * n + j * n + k] += Cs[i * n * n + k * n + j];
				Cs[i * n * n + k * n + j] = Cs[i * n * n + j * n + k];
			}
		}
		for (integer j = 0; j < n * n; j++)
		{
			Cs[i * n * n + j] *= 0.1;
		}
		for (integer j = 0; j < n; j++)
		{
			Cs[i * n * n + j * n + j] += n - j;
		}
	}

	CheckMemoryDeleted = new std::map<integer *, integer>;

	testStieSoftICA(Cs, n, p, N);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] Cs;

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

void testStieSoftICA(double *Cs, integer n, integer p, integer N, double *X, double *Xopt)
{
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	genrandseed(tt);

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	StieVariable StieX(n, p);
	StieX.RandInManifold();

	// Define the manifold
	Stiefel Domain(n, p);
	Domain.ChooseStieParamsSet2();
	Domain.SetHasHHR(false);

	// Define the Brockett problem
	StieSoftICA Prob(Cs, n, p, N);
	Prob.SetDomain(&Domain);

	Domain.CheckParams();

	//Domain.CheckIntrExtr(&StieX);
	//Domain.CheckRetraction(&StieX);
	//Domain.CheckDiffRetraction(&StieX, true);
	//Domain.CheckLockingCondition(&StieX);
	//Domain.CheckcoTangentVector(&StieX);
	//Domain.CheckIsometryofVectorTransport(&StieX);
	//Domain.CheckIsometryofInvVectorTransport(&StieX);
	//Domain.CheckVecTranComposeInverseVecTran(&StieX);
	//Domain.CheckTranHInvTran(&StieX);
	//Domain.CheckHaddScaledRank1OPE(&StieX);

	//// Check gradient and Hessian
	//Prob.CheckGradHessian(&StieX);

	// test RSD
	printf("********************************Check all line search algorithm in RSD*****************************************\n");
	for (integer i = 0; i < INPUTFUN; i++)
	{
		RSD *RSDsolver = new RSD(&Prob, &StieX);
		RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RSDsolver->Debug = FINALRESULT;
		RSDsolver->Max_Iteration = 20;
		RSDsolver->CheckParams();
		RSDsolver->Run();
		delete RSDsolver;
	}

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
	//printf("********************************Check all Formulas in RCG*************************************\n")
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
	//for (integer i = 0; i < INPUTFUN; i++)//LSALGOLENGTH
	//{
	//	RSD *RBFGSsolver = new RSD(&Prob, &StieX);
	//	RBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RBFGSsolver->Debug = ITERRESULT;//-- FINALRESULT;
	//	//RBFGSsolver->OutputGap = 20;
	//	//RBFGSsolver->LS_alpha = 0.1;
	//	//RBFGSsolver->InitSteptype = ONESTEP;
	//	RBFGSsolver->LS_ratio1 = 1.0/64;
	//	//RBFGSsolver->LS_ratio2 = 0.25;
	//	RBFGSsolver->CheckParams();
	//	RBFGSsolver->Run();

	//	delete RBFGSsolver;
	//}

	//RBFGS *RBFGSsolver = new RBFGS(&Prob, &StieX);
	//Domain.HasHHR = true;
	//RBFGSsolver->LineSearch_LS = WOLFE;
	//RBFGSsolver->Debug = ITERRESULT;
	//RBFGSsolver->OutputGap = 20;
	//RBFGSsolver->CheckParams();
	//RBFGSsolver->Run();
	//delete RBFGSsolver;

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

	//// test RTRNewton
	//printf("********************************Check RTRNewton*************************************\n");
	//RTRNewton RTRNewtonsolver(&Prob, &StieX);
	//RTRNewtonsolver.Debug = FINALRESULT;
	//RTRNewtonsolver.CheckParams();
	//RTRNewtonsolver.Run();
	//
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

	//// Check gradient and Hessian
	//Prob.CheckGradHessian(&StieX);
	//const Variable *xopt = RTRNewtonsolver.GetXopt();
	//Prob.CheckGradHessian(xopt);

	//const Variable *xopt = RTRNewtonsolver.GetXopt();
	//if (Xopt != nullptr)
	//{
	//	const double *xoptptr = xopt->ObtainReadData();
	//	for (integer i = 0; i < n * p; i++)
	//		Xopt[i] = xoptptr[i];
	//}
}

#endif

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs < 7)
    {
        mexErrMsgTxt("The number of arguments should be at least seven.\n");
    }
	double *Cs, *X, *Xopt;
	integer p, n, N, HasHHR, ParamSet;
	Cs = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	n = static_cast<integer> (mxGetScalar(prhs[2]));
	p = static_cast<integer> (mxGetScalar(prhs[3]));
	HasHHR = static_cast<integer> (mxGetScalar(prhs[4]));
	ParamSet = static_cast<integer> (mxGetScalar(prhs[5]));
	/* dimensions of input matrices */
	const mwSize *ptrdims = mxGetDimensions(prhs[0]);
	N = ptrdims[2];
    
    if(ptrdims[1] != n || ptrdims[0] != n)
    {
        mexErrMsgTxt("The size of matrix C is not correct.\n");
    }

	if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != p)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
    
	printf("(n, p, N):%d, %d, %d\n", n, p, N);

	///*create output matrix*/
	//plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL);
	//Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	//testStieSoftICA(Cs, n, p, N, X, Xopt);


	// Obtain an initial iterate by taking the Q factor of qr decomposition
	StieVariable *StieX = nullptr;

	StieX = new StieVariable(n, p);
	double *StieXptr = StieX->ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		StieXptr[i] = X[i];

	StieVariable *Stiesoln = nullptr;
	if (nrhs >= 8)
	{
		double *soln = mxGetPr(prhs[7]);
		Stiesoln = new StieVariable(n, p);
		double *Stiesolnptr = Stiesoln->ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
		{
			Stiesolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	Stiefel Domain(n, p);

	// Define the Brockett problem
	StieSoftICA Prob(Cs, n, p, N);
	Prob.SetDomain(&Domain);
	if (ParamSet == 1)
		Domain.ChooseStieParamsSet1();
	else
    if (ParamSet == 2)
		Domain.ChooseStieParamsSet2();
	else
    if (ParamSet == 3)
		Domain.ChooseStieParamsSet3();
	else
    if (ParamSet == 4)
		Domain.ChooseStieParamsSet4();
	else
    {
        printf("error: ParamSet parameter!");
        return;
    }
	Domain.SetHasHHR((HasHHR != 0));
	Domain.CheckParams();
	//Domain.CheckParams();
	ParseSolverParamsAndOptimizing(prhs[6], &Prob, StieX, Stiesoln, plhs);
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete StieX;
	delete Stiesoln;
	return;
}

#endif
