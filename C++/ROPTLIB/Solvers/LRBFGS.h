/*
This file defines the class of the limited-memory Riemannian BFGS method in [HGA2015]
	[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization. 
	SIAM Journal on Optimization, 25(3):1660?685, 2015

Solvers --> QuasiNewton --> SolversLS --> LRBFGS

---- WH
*/

#ifndef LRBFGS_H
#define LRBFGS_H

#include <cstring>
#include "Solvers/SolversLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRBFGS : public SolversLS{
	public:
		/*The contructor of LRBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research*/
		LRBFGS(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Destructor. Delete the arrays and vectors used in LRBFGS, i.e., series S and Y, and series RHO*/
		virtual ~LRBFGS();

		/*Check whether the parameters about LRBFGS are legal or not.*/
		virtual void CheckParams();

		/*Run the algorithm. New memory for S, Y and RHO. Then call SolversLS::Run*/
		virtual void Run();

		/*Call Solvers::SetProbX function and set up the temporary objects for LRBFGS algorithm.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

	protected:

		/*Compute the search direction. [HGA2015, Steps 3 to 14 in Algorithm 2]
			[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685, 2015.
			*/
		virtual void GetSearchDir();

		/*update the pairs of s and y. Add the latest one and remove the oldest one if necessary.
		transport them to the tangent space of x2*/
		virtual void UpdateData();

		/*Print information specific to LRBFGS*/
		virtual void PrintInfo();

	};
}; /*end of ROPTLIB namespace*/
#endif // end of RBROYDENFAMILY_H
