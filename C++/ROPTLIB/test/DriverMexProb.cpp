
#include "test/DriverMexProb.h"

#ifdef MATLAB_MEX_FILE

using namespace ROPTLIB;

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	CheckMemoryDeleted = new std::map<integer *, integer>;

	DriverMexProb(nlhs, plhs, nrhs, prhs);
	// check memory
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	return;
}

void DriverMexProb(int &nlhs, mxArray ** &plhs, int &nrhs, const mxArray ** &prhs)
{
	if (nrhs < 5)
	{
		mexErrMsgTxt("The number of arguments should be at least five.\n");
	}
	// Argument Checking:
	// First to third arguments should be function handles
	if (!mxIsClass(prhs[0], "function_handle") || !mxIsClass(prhs[1], "function_handle") || !mxIsClass(prhs[2], "function_handle"))
	{
		mexErrMsgTxt("At least one of first to third input arguments is not a function handle.");
	}
	// Fourth to fifth arguments are structures
	// The fourth one is SolverParams
	// The fifth one is ManiParams
	if (!mxIsStruct(prhs[3]) || !mxIsStruct(prhs[4]))
	{
		mexErrMsgTxt("At least one of fourth to fifth input arguments is not a structure.");
	}

	// Obtain manifold and iterate structure
	Manifold *domain = nullptr, **manifolds = nullptr;
	Variable *initialX = nullptr, *soln = nullptr;
	Element **elements;
	integer *powsinterval, numoftype, numoftotal;

	if (!ParseManiParams(prhs[4], manifolds, elements, powsinterval, numoftype, numoftotal))
	{
		mexErrMsgTxt("Parsing ManiParams fails.");
	}

	domain = new ProductManifold(manifolds, numoftype, powsinterval, numoftotal);
	mxArray *tmp = mexProblem::GetFieldbyName(prhs[4], 0, "IsCheckParams");
	if (tmp != nullptr)
	{
		if (fabs(mxGetScalar(tmp)) > std::numeric_limits<double>::epsilon()) // if the value is nonzero
		{
			domain->CheckParams();
		}
	}
	bool HasHHR = false;
	if (nrhs >= 6)
	{
		if (!mxIsDouble(prhs[5]))
		{
			mexErrMsgTxt("Sixth input argument is not a scalar.");
		}

		HasHHR = (static_cast<integer> (mxGetScalar(prhs[5])) != 0);
	}
	domain->SetHasHHR(HasHHR);
	initialX = new ProductElement(elements, numoftotal, powsinterval, numoftype);
	//initialX->Print("initialX");

	// initialize the initial iterate
	if (nrhs >= 7)
	{
		if (!mxIsStruct(prhs[6]))
		{
			mexErrMsgTxt("Seventh input argument is not a structure.");
		}
		mexProblem::ObtainElementFromMxArray(initialX, prhs[6]);
	}
	else
	{
		initialX->RandInManifold();
	}

	// 	initialX->Print("initialX", false);

	if (nrhs >= 8)
	{
		soln = new ProductElement(elements, numoftotal, powsinterval, numoftype);
		if (!mxIsStruct(prhs[7]))
		{
			mexErrMsgTxt("Seventh input argument is not a structure.");
		}
		mexProblem::ObtainElementFromMxArray(soln, prhs[7]);
	}
	//mexProblem::ObtainMxArrayFromElement(plhs[0], initialX);

	// Define the problem
	Problem *Prob = new mexProblem(prhs[0], prhs[1], prhs[2]);
	Prob->SetDomain(domain);
	//	Vector *egf = initialX->ConstructEmpty();
	//	Prob->EucGrad(initialX, egf);
	//	delete egf;
	// solve the optimization problem
	ParseSolverParamsAndOptimizing(prhs[3], Prob, initialX, soln, plhs);

	delete Prob;
	delete domain;
	delete initialX;
	delete soln;

	for (integer i = 0; i < numoftype; i++)
	{
		delete manifolds[i];
		delete elements[i];
	}
	delete[] manifolds;
	delete[] elements;
	delete[] powsinterval;
};

#endif // end of MATLAB_MEX_FILE
