/*
This file defines the Riemannian gradient sampling method, which is for optimizing partly smooth functions on Riemannian manifolds
See [Sec 7.2, Hua2014] or [HU2016] for details
[Hua2014]: Wen Huang, Optimization algorithm on Riemannian manifolds with applications, PhD thesis 2014.
[HU2016]: S. Hosseini, A. Uschmajew, A Riemannian gradient sampling algorithm for nonsmooth optimization on manifolds, INS Preprint No. 1607, 2016.

Solvers --> QuasiNewton --> SolversLS --> RGS

---- WH
*/

#ifndef RGS_H
#define RGS_H

#include "Others/MinPNormConHull.h"
#include "Solvers/SolversLSLPSub.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Compute min_{y in convex hull of gfs and prefgs are tangent vectors at the tangent space at x} ||y||
	It is defined in SphereConvexHull.h and SphereConvexHull.cpp */
	extern double MinNormConHull(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, Vector *Soln);

	class RGS : public SolversLSLPSub{
	public:
		/*The contructor of RGS method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		RGS is the initial inverse Hessian approximation. If the input is nullptr, then the identity is used as the initial approximation.
		insoln is the true solution. It is not required and only used for research.*/
		RGS(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*Destructor. Delete the vectors RGS*/
		virtual ~RGS();

	protected:

		/*Compute the search direction */
		void GetSearchDir(void);

	};
}; /*end of ROPTLIB namespace*/
#endif // end of RGS_H
