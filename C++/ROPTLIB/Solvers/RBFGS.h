/*
This file defines the class of the Riemannian BFGS method in [HGA2015]
[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
SIAM Journal on Optimization, 25(3):1660?685, 2015

Solvers --> QuasiNewton --> SolversLS --> RBFGS

---- WH
*/

#ifndef RBFGS_H
#define RBFGS_H

#include "Solvers/SolversLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class RBFGS : public SolversLS{
	public:
		/*The contructor of RBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		initialH is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		RBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH = nullptr, const Variable *insoln = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGS, i.e., s and y, H and tildeH*/
		virtual ~RBFGS();

		/*Check whether the parameters about RBFGS are legal or not.*/
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
		void GetSearchDir(void);

		/*Update the Hessian approximation if necessary*/
		void UpdateData(void);

		/*Print information specific to RBFGS*/
		virtual void PrintInfo();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of RBFGS_H
