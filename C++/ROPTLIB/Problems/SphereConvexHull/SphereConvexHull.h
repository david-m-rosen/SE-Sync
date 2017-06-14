/*
This file defines the class for 
min_{X \in S^{n-1}} tr((X.^2)^T W^T P W (X.^2)), where P is a N by N symmetric positive definite matrix
and W is a N by n matrix.

Problem --> SphereConvexHull

---- WH
*/

#ifndef SPHERECONVEXHULL_H
#define SPHERECONVEXHULL_H

#include "Manifolds/Sphere/Sphere.h"
#include "Manifolds/Sphere/SphereVariable.h"
#include "Manifolds/Sphere/SphereVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Solvers/SolversLS.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/LRBFGS.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereConvexHull : public Problem{
	public:
		/*W is an array of vectors, which are tangent vectors in the tangent space at x, lengthW is the length of W.
		HvRBFGSSub defines a linear mapping: P: v --> Pv. The manifold, Mani, of x and W is not necessary the same as the domain of this problem.*/
		SphereConvexHull(const Manifold *Mani, Variable *x, Vector **W, integer lengthW, QuasiNewton *solver, void (QuasiNewton::*HvRBFGSSub)(Vector *v, Vector *result));
		virtual ~SphereConvexHull();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		const Manifold *Mani;
		Variable *x;
		Vector **W;
		integer lengthW;
		QuasiNewton *solver;
		void (QuasiNewton::*Hv)(Vector *v, Vector *result);
	};

}; /*end of ROPTLIB namespace*/
#endif // end of GRASSRQ_H
