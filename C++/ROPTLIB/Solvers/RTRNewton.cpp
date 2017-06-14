
#include "Solvers/RTRNewton.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRNewton::RTRNewton(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void RTRNewton::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversTR::SetProbX(prob, initialx, insoln);
		prob->SetUseGrad(true);
		prob->SetUseHess(true);
	};

	void RTRNewton::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		SolverName.assign("RTRNewton");
	};

	void RTRNewton::HessianEta(Vector *Eta, Vector *result)
	{
		Prob->HessianEta(x1, Eta, result);
	};
}; /*end of ROPTLIB namespace*/
