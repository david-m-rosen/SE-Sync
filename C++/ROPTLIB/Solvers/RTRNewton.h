/*
This file defines the class of the Riemannian trust-region Newton method in [ABG2007]
	[ABG2007]: P.-A. Absil, C. G. Baker, and K. A. Gallivan. Trust-region methods on Riemannian manifolds. 
		Foundations of Computational Mathematics, 7(3):303?30, 2007.

Solvers --> SolversTR --> RTRNewton

---- WH
*/

#ifndef RTRNEWTON_H
#define RTRNEWTON_H

#include "Solvers/SolversTR.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RTRNewton : public SolversTR{
	public:
		/*The contructor of RTRNewton method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		RTRNewton(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Call Solvers::SetProbX function; initialize temporary vectors; and indicate RTRNewton needs action of Hessian.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();
	protected:
		/*Compute result = H[Eta], where H is the Hessian*/
		virtual void HessianEta(Vector *Eta, Vector *result);
	};
}; /*end of ROPTLIB namespace*/

#endif // end of RTRNEWTON_H