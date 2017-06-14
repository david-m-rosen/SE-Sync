
#include "test/TestStieBrockett.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSTIEBROCKETT is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTSTIEBROCKETT)

/*Help to check the memory leakage problem. No necesary any more.*/
std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	/*Set the random seed*/
	unsigned tt = (unsigned)time(NULL);
	tt = 0;
	genrandseed(tt);

	// size of the Stiefel manifold
	integer n = 128, p = 8;

	// Generate the matrices in the Brockett problem.
	double *B = new double[n * n + p];
	double *D = B + n * n;
	/*B is an n by n matrix*/
	for (integer i = 0; i < n; i++)
	{
		for (integer j = i; j < n; j++)
		{
			B[i + j * n] = genrandnormal();
			B[j + i * n] = B[i + j * n];
		}
	}
	/*D is a diagonal matrix.*/
	for (integer i = 0; i < p; i++)
		D[i] = static_cast<double> (1);//-- (i + 1);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	testStieBrockett(B, D, n, p);
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

	const StieBrockett *P = dynamic_cast<StieBrockett *> (const_cast<Problem*> (prob));
	const double *etaTV = eta1->ObtainReadData();
	const double *xM = x1->ObtainReadData();

	integer n = P->n, p = P->p, length = n * p;
	double denor, nume;
	double *B = P->B, *D = P->D;
	double *tmp = new double [n * p];
	dgemm_(GLOBAL::N, GLOBAL::N, &n, &p, &n, &GLOBAL::DONE, B, &n, const_cast<double *> (etaTV), &n, &GLOBAL::DZERO, tmp, &n);

	denor = ddot_(&length, tmp, &GLOBAL::IONE, const_cast<double *> (etaTV), &GLOBAL::IONE);
	nume = ddot_(&length, tmp, &GLOBAL::IONE, const_cast<double *> (xM), &GLOBAL::IONE);

	delete[] tmp;
	return (- nume / denor < 0) ? 1 : - nume / denor;
}

void testStieBrockett(double *B, double *D, integer n, integer p, double *X, double *Xopt)
{
	//// choose a random seed
	//unsigned tt = (unsigned)time(NULL);
	////tt = 0;
	//genrandseed(tt);
	StieVariable StieX(n, p);

	if (X == nullptr)
	{/*If X is not defined before, then obtain an initial iterate by taking the Q factor of qr decomposition*/
		StieX.RandInManifold();
	}
	else
	{/*Otherwise, using the input orthonormal matrix as the initial iterate*/
		double *StieXptr = StieX.ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
			StieXptr[i] = X[i];
	}

	// Define the manifold
	Stiefel Domain(n, p);
	//Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);
	/*The domain of the problem is a Stiefel manifold*/
	Prob.SetDomain(&Domain);

	/*Output the parameters of the domain manifold*/
	Domain.CheckParams();

	//Domain.CheckRetraction(&StieX);
	//Domain.CheckDiffRetraction(&StieX);
	//Domain.CheckLockingCondition(&StieX);
	//Domain.CheckcoTangentVector(&StieX);
	//Domain.CheckIsometryofVectorTransport(&StieX);
	//Domain.CheckIsometryofInvVectorTransport(&StieX);
	//Domain.CheckVecTranComposeInverseVecTran(&StieX);
	//Domain.CheckTranHInvTran(&StieX);
	//Domain.CheckHaddScaledRank1OPE(&StieX);

	//RBFGSLPSub *RBFGSLPSubsolver = new RBFGSLPSub(&Prob, &StieX);
	//RBFGSLPSubsolver->Debug = FINALRESULT;
	//RBFGSLPSubsolver->OutputGap = 10;
	//RBFGSLPSubsolver->lambdaLower = 1e-3;
	//RBFGSLPSubsolver->lambdaUpper = 1e3;
	//RBFGSLPSubsolver->CheckParams();
	//RBFGSLPSubsolver->Run();
	//delete RBFGSLPSubsolver;
	//return;

	//RBFGSLPSubsolver = new RBFGSLPSub(&Prob, &StieX);
	//RBFGSLPSubsolver->Debug = FINALRESULT;
	//RBFGSLPSubsolver->OutputGap = 10;
	//RBFGSLPSubsolver->lambdaLower = 1e-7;
	//RBFGSLPSubsolver->lambdaUpper = 1e7;
	//RBFGSLPSubsolver->CheckParams();
	//RBFGSLPSubsolver->Run();
	//delete RBFGSLPSubsolver;


	//RSD *RSDsolver = new RSD(&Prob, &StieX);
	//RSDsolver->Debug = FINALRESULT;
	//RSDsolver->InitSteptype = ONESTEP;
	//RSDsolver->Max_Iteration = 2000;
	//RSDsolver->CheckParams();
	//RSDsolver->Run();
	//delete RSDsolver;

	//RBFGS *RBFGSsolver = new RBFGS(&Prob, &StieX);
	//RBFGSsolver->Debug = FINALRESULT;
	//RBFGSsolver->CheckParams();
	//RBFGSsolver->Run();
	//delete RBFGSsolver;

	//// test RSD
	//printf("********************************Check all line search algorithm in RSD*****************************************\n");
	//for (integer i = 0; i < INPUTFUN; i++)
	//{
	//	RSD *RSDsolver = new RSD(&Prob, &StieX);
	//	RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RSDsolver->Debug = FINALRESULT;
	//	RSDsolver->Max_Iteration = 2000;
	//	RSDsolver->CheckParams();
	//	RSDsolver->Run();
	//	delete RSDsolver;
	//}

	//// test RNewton
	//printf("********************************Check all line search algorithm in RNewton*************************************\n");
	//for (integer i = INPUTFUN; i < LSALGOLENGTH; i++)
	//{
	//	RNewton *RNewtonsolver = new RNewton(&Prob, &StieX);
	//	RNewtonsolver->LineSearch_LS = static_cast<LSAlgo> (i);
	//	RNewtonsolver->Debug = ITERRESULT;
	//	/*Uncomment following two lines to use the linesearch algorithm defined by the function "LinesearchInput".*/
	//	RNewtonsolver->LineSearch_LS = INPUTFUN;
	//	RNewtonsolver->LinesearchInput = &LinesearchInput;
	//	RNewtonsolver->Max_Iteration = 100;
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
	//	RCGsolver->LineSearch_LS = ARMIJO;
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

	// test LRBFGS
	printf("********************************Check all line search algorithm in LRBFGS*************************************\n");
	//Domain.ChooseStieParamsSet4();
	for (integer i = 0; i < 1; i++)//INPUTFUN
	{
		LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, &StieX);
		LRBFGSsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		LRBFGSsolver->Debug = FINALRESULT; //ITERRESULT;// 
		LRBFGSsolver->OutputGap = 1;
		LRBFGSsolver->Max_Iteration = 1000;
		LRBFGSsolver->Accuracy = 1e-6;
		LRBFGSsolver->Tolerance = 1e-6;
		LRBFGSsolver->Finalstepsize = 1;
		LRBFGSsolver->Num_pre_funs = 0;
		LRBFGSsolver->BBratio = 1;
		LRBFGSsolver->Num_pre_BB = 0;
		LRBFGSsolver->InitSteptype = ONESTEP;
		LRBFGSsolver->LengthSY = 0;
		LRBFGSsolver->CheckParams();
		LRBFGSsolver->Run();
		//// Check the correctness of gradient and Hessian at the initial iterate
		//Prob.CheckGradHessian(&StieX);
		//const Variable *xopt = LRBFGSsolver->GetXopt();
		//// Check the correctness of gradient and Hessian at the final iterate of RTRNewton method
		//Prob.CheckGradHessian(xopt);
		delete LRBFGSsolver;
	}

	//// test RTRSD
	//printf("********************************Check RTRSD*************************************\n");
	//RTRSD RTRSDsolver(&Prob, &StieX);
	//RTRSDsolver.Debug = FINALRESULT;
	//RTRSDsolver.Max_Iteration = 5000;
	//RTRSDsolver.CheckParams();
	//RTRSDsolver.Run();

	//// test RTRNewton
	//printf("********************************Check RTRNewton*************************************\n");
	//RTRNewton RTRNewtonsolver(&Prob, &StieX);
	//RTRNewtonsolver.Debug = FINALRESULT;
	//RTRNewtonsolver.CheckParams();
	//RTRNewtonsolver.Run();

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

	//// Check the correctness of gradient and Hessian at the initial iterate
	//Prob.CheckGradHessian(&StieX);
	//const Variable *xopt = RTRNewtonsolver.GetXopt();
	//// Check the correctness of gradient and Hessian at the final iterate of RTRNewton method
	//Prob.CheckGradHessian(xopt);
    
	////Output the optimizer obtained by RTRNewton method
	//if (Xopt != nullptr)
	//{
	//	const double *xoptptr = xopt->ObtainReadData();
	//	for (integer i = 0; i < n * p; i++)
	//		Xopt[i] = xoptptr[i];
	//}
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
	if(nrhs < 6)
	{
		mexErrMsgTxt("The number of arguments should be at least six.\n");
	}
	double *B, *D, *X, *Xopt;
	B = mxGetPr(prhs[0]);
	D = mxGetPr(prhs[1]);
	X = mxGetPr(prhs[2]);
	/* dimensions of input matrices */
	integer p, n, HasHHR, Paramset;
	n = mxGetM(prhs[0]);
	p = mxGetM(prhs[1]);

	/*Check the correctness of the inputs*/
	if(mxGetN(prhs[0]) != n)
	{
		mexErrMsgTxt("The size of matrix B is not correct.\n");
	}
	if(mxGetN(prhs[1]) != 1)
	{
		mexErrMsgTxt("The size of the D is not correct!\n");
	}
	if(mxGetM(prhs[2]) != n || mxGetN(prhs[2]) != p)
	{
		mexErrMsgTxt("The size of the initial X is not correct!\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[3]));
	Paramset = static_cast<integer> (mxGetScalar(prhs[4]));

	printf("(n, p):%d,%d\n", n, p);

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n, p, mxREAL);
	Xopt = mxGetPr(plhs[0]);

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;
	//	testStieBrockett(B, D, n, p, X, Xopt);

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	StieVariable StieX(n, p);
	double *StieXptr = StieX.ObtainWriteEntireData();
	for (integer i = 0; i < n * p; i++)
		StieXptr[i] = X[i];

	StieVariable *Stiesoln = nullptr;
	if (nrhs >= 7)
	{
		double *soln = mxGetPr(prhs[6]);
		Stiesoln = new StieVariable(n, p);
		double *Stiesolnptr = Stiesoln->ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
		{
			Stiesolnptr[i] = soln[i];
		}
	}
	// Define the manifold
	Stiefel Domain(n, p);
	if (Paramset == 1)
		Domain.ChooseStieParamsSet1();
	else if (Paramset == 2)
		Domain.ChooseStieParamsSet2();
	else if (Paramset == 3)
		Domain.ChooseStieParamsSet3();

	// Define the Brockett problem
	StieBrockett Prob(B, D, n, p);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[5], &Prob, &StieX, Stiesoln, plhs);

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete Stiesoln;
	return;
}

#endif
