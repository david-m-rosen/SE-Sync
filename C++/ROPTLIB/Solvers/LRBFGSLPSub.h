/*
This file defines the class of the line search algorithm for locally lipschitz functions on Riemannian manifolds

Solvers --> QuasiNewton --> SolversLS --> SolversLPSub --> LRBFGSLPSub

---- WH
*/

#ifndef LRBFGSLPSUB_H
#define LRBFGSLPSUB_H

#include "Solvers/SolversLS.h"
#include "Problems/SphereConvexHull/SphereConvexHull.h"
#include "Manifolds/Sphere/Sphere.h"
#include "Solvers/LRBFGS.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/SolversLSLPSub.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRBFGSLPSub : public SolversLSLPSub{
	public:
		/*Run the algorithm. This function gives the framework for the linesearch method*/
		virtual void Run(void);

		/*The contructor of RBFGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research*/
		LRBFGSLPSub(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGSLPSub, i.e., s and y, H and tildeH*/
		virtual ~LRBFGSLPSub();

		/*Check whether the parameters about RBFGSLPSub are legal or not.*/
		virtual void CheckParams();

		/*Initialize the solvers by calling the "SetProbX" and "SetDefultParams" functions.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void Initialization(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2, s, y, H and tildeH and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

	protected:

		/*Print information specific to SolversLPSub*/
		virtual void PrintInfo();

		/*Update the Hessian approximation if necessary*/
		void UpdateData(void);
	};
}; /*end of ROPTLIB namespace*/
#endif // end of RBFGSLPSUB_H
