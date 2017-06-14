/*
This file defines the abstract base class for problem classes. All the problem classes
should be derived from this class.

Problem

---- WH
*/

#ifndef PROBLEM_H
#define PROBLEM_H

//#include <cmath>

#include "Manifolds/Element.h"
#include "Manifolds/Manifold.h"
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Declaration of Manifold*/
	class Manifold;

	class Problem{
	public:
		/*This function indicates this is an abstract class.*/
		virtual ~Problem(void) = 0;

		/*Evaluate the cost function at iterate x. It must be overloaded by derived class*/
		virtual double f(Variable *x) const = 0;

		/* It calls "RieGrad" and obtain the extrinsic representation of the Riemannian gradient.
		After that this function may or may not convert the extrinsic representation to intrinsic representation
		based on the "IsIntrApproach" in the domain manifold. */
		virtual void Grad(Variable *x, Vector *gf) const;

		/*Compute the action of the Riemannian Hessian of the cost function at iterate x, i.e., xix = Hess f(x) [etax].
		It calls "RieHess" and obtain the extrinsic representation of the Action of the Hessian.
		After that this function may or may not convert the extrinsic representation to intrinsic representation
		based on the "IsIntrApproach" in the domain manifold. */
		virtual void HessianEta(Variable *x, Vector *etax, Vector *xix) const;

		/*Compute the Riemanian  gradient of the cost function at iterate x.
		User can override this function. Otherwise, this function will call "EucGrad" to
		obtain the Euclidean gradient. Then convert the Euclidean gradient to Riemannian gradient.*/
		virtual void RieGrad(Variable *x, Vector *gf) const;

		/*Compute the action of the Riemanian Hessian of the cost function at iterate x.
		etax and xix must be different arguments, i.e.,
		calling this function by
		RieHessianEta(x, xix, xix);
		is illegal.
		User can override this function. Otherwise, this function will call "EucHessianEta" to
		obtain the Euclidean action of the Hessian. Then convert the Euclidean action of the Hessian
		to Riemannian action of the Hessian. */
		virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		/*Compute the Euclidean gradient of f*/
		virtual void EucGrad(Variable *x, Vector *egf) const;

		/*Compute the action of the Euclidean Hessian,
		etax and xix must be different arguments, i.e.,
		calling this function by
		EucHessianEta(x, xix, xix);
		is illegal. */
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		/*Check the correctness of the Riemannian gradient and Riemannian Hessian.
		See details in the user manual*/
		virtual void CheckGradHessian(const Variable *x) const;

		/*Set the domain of the cost function.*/
		virtual void SetDomain(Manifold *inDomain);

		/*Obtain the domain manifold of the cost function*/
		inline Manifold *GetDomain(void) const { return Domain; };

		/*Mark whether the gradient is used in the algorithm or not.*/
		virtual void SetUseGrad(bool usegrad) const;

		/*Mark whether the action of the Hessian is used in the algorithm or not.*/
		virtual void SetUseHess(bool usehess) const;

		/*Get whether the gradient is used*/
		inline bool GetUseGrad(void) const { return UseGrad; };

		/*Get whether the action of the Hessian is used*/
		inline bool GetUseHess(void) const { return UseHess; };
	protected:
		Manifold *Domain; /*Domain of the cost function. It is required to be assigned.*/
		mutable bool UseGrad; /*Mark whether the gradient is used*/
		mutable bool UseHess; /*Mark whether the action of the Hessian is used.*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of PROBLEM_H
