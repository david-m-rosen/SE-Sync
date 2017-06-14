
#include "Solvers/MRankAdaptive.h"

/*Define the namespace*/
namespace ROPTLIB{

	MRankAdaptive::MRankAdaptive(const Problem *prob, const Variable *initialx, Solvers *inSolver, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
		innerSolver = inSolver;
	};

	void MRankAdaptive::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversLS::SetProbX(prob, initialx, insoln);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void MRankAdaptive::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		SolverName.assign("MRankAdaptive");
		// setup the default values for parameters
		//Reps1 = ;
		//Reps2 = ;
		//Reps3 = ;
		//Rca = ;
		//Rcr = ;
		//Rtau2 = ;
		//RDelta = ;
	};

	void MRankAdaptive::GetSearchDir(void)
	{
		Mani->ScaleTimesVector(x1, -1, gf1, eta1);
	};

	void MRankAdaptive::InitialStepSize(void)
	{
		if (iter == 0)
			stepsize = Initstepsize / ngf;
		else
		{
			stepsize = 1.01 * 2.0 * (f1 - pre_funs.front()) / initialslope;
			stepsize = (stepsize < std::numeric_limits<double>::epsilon()) ? Initstepsize / ngf : stepsize;
		}
	};

	void MRankAdaptive::Run(void)
	{
		// implement the main algorithm
		// 1), the line search method can use the function in SolversLS
		// 2), the inner solver is the innerSolver, which can be any Riemannian solvers or other solvers that is derived from Solvers.h
		// 3), The stopping criterion of innerSolver is set up by using InnerStop function:
		innerSolver->StopPtr = &MRA::InnerStop;
		// 4), The Rank-related-retraction need to be implemented in the Lowrank.h
		// 5), Since the rank is variant during iteration, please make sure the manifold, variable and vector need to be have consistant size.
	};

	namespace MRA{
		bool InnerStop(Variable *x, Vector *gf, double f, double ngf, double ngf0, const Problem *prob, const Solvers *solver)
		{
			// If either ngf is small enough using Reps3 or the singular values of x has significant bias using MRA::RDelta,
			// then stop.
			return true;
		};
	};
}; /*end of ROPTLIB namespace*/
