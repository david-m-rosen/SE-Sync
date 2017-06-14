/*
This file defines the class of the Riemannian trust-region Newton method in [ABG2007]
[ABG2007]: P.-A. Absil, C. G. Baker, and K. A. Gallivan. Trust-region methods on Riemannian manifolds.
Foundations of Computational Mathematics, 7(3):303?30, 2007.

Solvers --> SolversTR --> RTRSD

---- WH
*/

#ifndef RTRSD_H
#define RTRSD_H

#include "Solvers/SolversTR.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RTRSD : public SolversTR{
	public:
		/*The contructor of RTRSD method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		RTRSD(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Call Solvers::SetProbX function; initialize temporary vectors; and indicate RTRSD does not need action of Hessian.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();
	protected:
		/*Set result = Eta. In other words, the Hessian approximation is just an identity*/
		virtual void HessianEta(Vector *Eta, Vector *result);
	};
}; /*end of ROPTLIB namespace*/

#endif // end of RTRSD_H
