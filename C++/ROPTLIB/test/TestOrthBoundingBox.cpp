
#include "test/TestOrthBoundingBox.h"

using namespace ROPTLIB;

/*If the file is not compiled in Matlab and TESTSTIEBROCKETT is defined in def.h file, then using the following
main() function as the entrance. */
#if !defined(MATLAB_MEX_FILE) && defined(TESTORTHBOUNDINGBOX)

/*Help to check the memory leakage problem. No necesary any more.*/
std::map<integer *, integer> *CheckMemoryDeleted;

int main(void)
{
	/*Set the random seed*/
	unsigned tt = (unsigned)time(NULL);
	printf("random seed:%ud\n", tt);
	tt = 0;
	//tt = 1463379963;
	//tt = 1462047922;
	genrandseed(tt);

	// size of the problem
	integer d = 3, n = 10;

	// Generate the matrices in the Brockett problem.
	double *E = new double[d * n];
	/*E is an d by n matrix*/
	//for (integer i = 0; i < n * d; i++)
	//{
	//	E[i] = genrandreal();
	//	//E[i] = genrandnormal();
	//}

	CheckMemoryDeleted = new std::map<integer *, integer>;
	for (integer i = 3; i < 4; i++)
	{
		//genrandseed(i);
		testOrthBoundingBox(E, d, n);
	}
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete[] E;
#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
}

#endif

/*The main test function*/
void testOrthBoundingBox(double *E, integer d, integer n, double *X, double *Xopt)
{
	//// choose a random seed
	//unsigned tt = (unsigned)time(NULL);
	////tt = 0;
	//genrandseed(tt);
	for (integer i = 0; i < n * d; i++)
	{
		E[i] = genrandreal();
	}
	OrthGroupVariable OrthX(d);

	if (X == nullptr)
	{/*If X is not defined before, then obtain an initial iterate by taking the Q factor of qr decomposition*/
		OrthX.RandInManifold();
	}
	else
	{/*Otherwise, using the input orthonormal matrix as the initial iterate*/
		double *OrthXptr = OrthX.ObtainWriteEntireData();
		for (integer i = 0; i < d * d; i++)
			OrthXptr[i] = X[i];
	}

	// Define the manifold
	OrthGroup Domain(d);
	Domain.SetHasHHR(true); /*set whether the manifold uses the idea in [HGA2015, Section 4.3] or not*/

	// Define the Bounding box problem
	OrthBoundingBox Prob(E, d, n);
	/*The domain of the problem is an orthogonal group manifold*/
	Prob.SetDomain(&Domain);

	/*Output the parameters of the domain manifold*/
	//Domain.CheckParams();

	//Domain.CheckRetraction(&OrthX);
	//Domain.CheckDiffRetraction(&OrthX);
	//Domain.CheckLockingCondition(&OrthX);
	//Domain.CheckcoTangentVector(&OrthX);
	//Domain.CheckIsometryofVectorTransport(&OrthX);
	//Domain.CheckIsometryofInvVectorTransport(&OrthX);
	//Domain.CheckVecTranComposeInverseVecTran(&OrthX);
	//Domain.CheckTranHInvTran(&OrthX);
	//Domain.CheckHaddScaledRank1OPE(&OrthX);

	//Prob.CheckGradHessian(&OrthX);

	printf("RBFGSSub\n");
	RBFGSLPSub * RBFGSLPSubsolver = new RBFGSLPSub(&Prob, &OrthX);
	RBFGSLPSubsolver->Debug = ITERRESULT;
	//RBFGSLPSubsolver->Tolerance = 1e-12;
	//RBFGSLPSubsolver->InitSteptype = QUADINTMOD;
	//RBFGSLPSubsolver->Maxstepsize = 1e6;
	RBFGSLPSubsolver->NumExtraGF = Domain.GetIntrDim() * 10;
	RBFGSLPSubsolver->lambdaLower = 1e-16;
	RBFGSLPSubsolver->lambdaUpper = 1e16;
	RBFGSLPSubsolver->Max_Iteration = 100;
	RBFGSLPSubsolver->OutputGap = 1;
	RBFGSLPSubsolver->CheckParams();
	RBFGSLPSubsolver->Run();

	//printf("LRBFGSSub\n");
	//LRBFGSLPSub *LRBFGSLPSubsolver = new LRBFGSLPSub(&Prob, &OrthX, RBFGSLPSubsolver->GetXopt());
	//LRBFGSLPSubsolver->Debug = ITERRESULT;
	////LRBFGSLPSubsolver->Tolerance = 1e-12;
	////LRBFGSLPSubsolver->InitSteptype = QUADINTMOD;
	//LRBFGSLPSubsolver->NumExtraGF = 5;
	//LRBFGSLPSubsolver->lambdaLower = 1e-3;
	//LRBFGSLPSubsolver->lambdaUpper = 1e3;
	//LRBFGSLPSubsolver->Max_Iteration = 10000;
	//LRBFGSLPSubsolver->OutputGap = 100;
	////LRBFGSLPSubsolver->CheckParams();
	//LRBFGSLPSubsolver->Run();
	//for (integer i = 0; i < LRBFGSLPSubsolver->GetlengthSeries(); i++)
	//{
	//	std::cout << i << ":" << LRBFGSLPSubsolver->GetdistSeries()[i] << std::endl;
	//}
	//delete LRBFGSLPSubsolver;
	delete RBFGSLPSubsolver;

	//printf("RBFGS\n");
	//RBFGS *RBFGSsolver = new RBFGS(&Prob, &OrthX);
	//RBFGSsolver->LineSearch_LS = WOLFELP;
	//RBFGSsolver->Stop_Criterion = PSSUBGRAD;
	//RBFGSsolver->Debug = ITERRESULT;
	//RBFGSsolver->Max_Iteration = 1000;
	//RBFGSsolver->NumExtraGF = 5;
	//RBFGSsolver->OutputGap = 100;
	//RBFGSsolver->Minstepsize = 1e-7;
	////RBFGSsolver->CheckParams();
	//RBFGSsolver->Run();

	////Vector *subgf = Domain.GetEMPTYINTR()->ConstructEmpty();
	//////Prob.f(&OrthX);
	//////Prob.Grad(&OrthX, subgf);
	////printf("test\n");//----
	////Prob.f(const_cast<Variable *>(RBFGSsolver->GetXopt()));
	////Prob.Grad(const_cast<Variable *>(RBFGSsolver->GetXopt()), subgf);
	////subgf->Print("subgf:");//----

	////integer numdim = Domain.GetIntrDim() * 100;
	////Vector **dgfs = new Vector*[numdim];
	////for (integer i = 0; i < numdim; i++)
	////	dgfs[i] = Domain.GetEMPTYINTR()->ConstructEmpty();
	////integer length = 0;
	////Prob.AllDirectionDerivative(RBFGSsolver->GetXopt(), dgfs, numdim, 1e-6, length);
	//////Prob.AllDirectionDerivative(&OrthX, dgfs, numdim, 1e-4, length);

	////printf("length:" << length << std::endl;//---
	//////for (integer i = 0; i < length; i++)
	//////{
	//////	printf("i:" << i << std::endl;
	//////	dgfs[i]->Print("dgfs");//---
	//////}
	////printf("norm:" << MinPNormConHull(&Domain, const_cast<Variable *>(RBFGSsolver->GetXopt()), dgfs, length, nullptr, nullptr, 0) << std::endl;

	////for (integer i = 0; i < numdim; i++)
	////	delete dgfs[i];
	////delete[] dgfs;
	////delete subgf;
	//delete RBFGSsolver;

	//printf("RGS\n");
	//RGS *RGSsolver = new RGS(&Prob, &OrthX);
	//RGSsolver->LineSearch_LS = ARMIJO;
	//RGSsolver->NumExtraGF = 5;
	////RGSsolver->Eps = 1e-5;
	//RGSsolver->Max_Iteration = 500;
	//RGSsolver->Debug = ITERRESULT;
	//RGSsolver->OutputGap = 100;
	//RGSsolver->Run();
	//delete RGSsolver;

	////Output the optimizer obtained by RBFGSLPSub method
	//if (Xopt != nullptr)
	//{
	//	const double *xoptptr = xopt->ObtainReadData();
	//	for (integer i = 0; i < n * p; i++)
	//		Xopt[i] = xoptptr[i];
	//}
}

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
	if(nrhs < 5)
	{
		mexErrMsgTxt("The number of arguments should be at least five.\n");
	}
	double *E, *X;
	E = mxGetPr(prhs[0]);
	X = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer d, n, HasHHR, Paramset;
	d = mxGetM(prhs[0]);
	n = mxGetN(prhs[0]);

	/*Check the correctness of the inputs*/
	if(mxGetN(prhs[1]) != d || mxGetM(prhs[1]) != d)
	{
		mexErrMsgTxt("The sizes of matrix E and X do not match.\n");
	}
	HasHHR = static_cast<integer> (mxGetScalar(prhs[2]));
	Paramset = static_cast<integer> (mxGetScalar(prhs[3]));

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	// Obtain an initial iterate by taking the Q factor of qr decomposition
	OrthGroupVariable OrthX(d);
	double *OrthXptr = OrthX.ObtainWriteEntireData();
	for (integer i = 0; i < d * d; i++)
		OrthXptr[i] = X[i];

	OrthGroupVariable *Orthsoln = nullptr;
	if (nrhs >= 6)
	{
		double *soln = mxGetPr(prhs[5]);
		Orthsoln = new OrthGroupVariable(d);
		double *Orthsolnptr = Orthsoln->ObtainWriteEntireData();
		for (integer i = 0; i < d * d; i++)
		{
			Orthsolnptr[i] = soln[i];
		}
	}

	// Define the manifold
	OrthGroup Domain(d);
	if (Paramset == 1)
		Domain.ChooseStieParamsSet1();
	else if (Paramset == 2)
		Domain.ChooseStieParamsSet2();
	else if (Paramset == 3)
		Domain.ChooseStieParamsSet3();

	// Define the Brockett problem
	OrthBoundingBox Prob(E, d, n);
	Prob.SetDomain(&Domain);

	Domain.SetHasHHR(HasHHR != 0);
	//Domain.CheckParams();

	// Call the function defined in DriverMexProb.h
	ParseSolverParamsAndOptimizing(prhs[4], &Prob, &OrthX, Orthsoln, plhs);

	mexProblem::ObtainElementFromMxArray(&OrthX, plhs[0]);

	integer numdim = Domain.GetIntrDim() * 10;
	Vector **dgfs = new Vector*[numdim];
	for (integer i = 0; i < numdim; i++)
		dgfs[i] = Domain.GetEMPTYINTR()->ConstructEmpty();

	integer length = 0;
	double testzero = 0;

	/*Get all possible directional derivatives*/
	//Prob.AllDirectionDerivative(&OrthX, dgfs, numdim, 1e-6, length);
	//
	///*gradient sampling*/
	//integer length = numdim;
	//Domain.RandomTangentVectors(&OrthX, numdim, dgfs);
	//double tmp = 0, Eps = 1e-6;

	//for (integer i = 1; i < numdim; i++)
	//{
	//	tmp = sqrt(Domain.Metric(&OrthX, dgfs[i], dgfs[i]));
	//	Domain.ScaleTimesVector(&OrthX, genrandreal() * Eps / tmp, dgfs[i], dgfs[i]);
	//}
	//OrthGroupVariable OrthY(d);
	//Prob.f(&OrthX);
	//Prob.Grad(&OrthX, dgfs[0]);
	//for (integer i = 1; i < numdim; i++)
	//{
	//	Domain.Retraction(&OrthX, dgfs[i], &OrthY);
	//	Prob.f(&OrthY);
	//	Prob.Grad(&OrthY, dgfs[i]);
	//	//Domain.InverseVectorTransport(&OrthX, dgfs[i], &OrthY, dgfs[i], dgfs[i]);
	//}

	//testzero = MinPNormConHullRMethod(&Domain, &OrthX, dgfs, length, nullptr, nullptr, nullptr);
	//printf("length:" << length << ", norm:" << testzero << std::endl;

	plhs[15] = mxCreateDoubleScalar(testzero);
	for (integer i = 0; i < numdim; i++)
		delete dgfs[i];
	delete[] dgfs;

	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	delete Orthsoln;
	return;
}

#endif
