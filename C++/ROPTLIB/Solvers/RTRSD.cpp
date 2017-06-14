
#include "Solvers/RTRSD.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRSD::RTRSD(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void RTRSD::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversTR::SetProbX(prob, initialx, insoln);
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RTRSD::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		theta = 0.1;
		kappa = 0.9;
		SolverName.assign("RTRSD");
	};

	void RTRSD::HessianEta(Vector *Eta, Vector *result)
	{
		Eta->CopyTo(result);
	};
}; /*end of ROPTLIB namespace*/
