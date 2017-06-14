/*
This file defines the class of the Riemannian BFGS method in [RW0212]
	[RW0212]: W. Ring and B. Wirth. Optimization methods on Riemannian manifolds and their application to shape space. 
	SIAM Journal on Optimization, 22(2):596?27, January 2012.

Solvers --> QuasiNewton --> SolversLS --> RWRBFGS

---- WH
*/
#ifndef RWRBFGS_H
#define RWRBFGS_H

#include "Solvers/SolversLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RWRBFGS : public SolversLS{
	public:
		/*The contructor of RBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		RWRBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr, const Variable *insoln = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGS, i.e., s and y, H and tildeH*/
		virtual ~RWRBFGS();

		/*Check whether the parameters about RWRBFGS are legal or not.*/
		virtual void CheckParams();

		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr, const Variable *insoln = nullptr);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr, const Variable *insoln = nullptr);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();
	protected:
		/*Compute the search direction. eta1 = H (-gf1) */
		virtual void GetSearchDir();

		/*Update the Hessian approximation if necessary*/
		virtual void UpdateData();

		/*Print information specific to RWRBFGS*/
		virtual void PrintInfo();
	};
}; /*end of ROPTLIB namespace*/

#endif // end of RWRBFGS_H
