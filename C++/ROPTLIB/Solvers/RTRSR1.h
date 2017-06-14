/*
This file defines the class of the Riemannian trust-region Newton method in [HAG2015]
	[HAG2015]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method. 
	Mathematical Programming, 150(2):179?16, 2015

Solvers --> SolversTR --> RTRSR1

---- WH
*/

#ifndef RTRSR1_H
#define RTRSR1_H

#include "Solvers/SolversTR.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RTRSR1 : public SolversTR{
	public:
		/*The contructor of RTRSR1 method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialB is the initial Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		RTRSR1(const Problem *prob, const Variable *initialx, LinearOPE *initialB = nullptr, const Variable *insoln = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGS, i.e., s and y, H and tildeH*/
		virtual ~RTRSR1();

		/*Check whether the parameters about RBFGS are legal or not.*/
		virtual void CheckParams();

		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialB = nullptr, const Variable *insoln = nullptr);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialB is the initial Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialB = nullptr, const Variable *insoln = nullptr);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

	protected:
		/*Compute the action of the Hessian approximation. result = B [Eta]*/
		virtual void HessianEta(Vector *Eta, Vector *result);

		/*Update the Hessian approximation if necessary*/
		virtual void UpdateData();

		/*This function is called when the candidate is accepted. It transports the Hessian approximation B
		from the tangent space of x1 to x2*/
		virtual void Acceptence();

		/*Print information specific to RTRSR1*/
		virtual void PrintInfo();
	};
}; /*end of ROPTLIB namespace*/

#endif // end of RTRSR1_H
