/*
This file defines the class for the sphere in L^2([0, 1], R), a point x in L^2([0, 1], R) is represented by its values on uniformly-spaced
grid, i.e., ( x(0/(n-1)), x(1/(n-1)), cdots, x((n-1)/(n-1)) ). Composite trapezoidal rule is used to define its Riemannian metric.

Manifold --> L2Sphere

---- WH
*/

#ifndef L2SPHERE_H
#define L2SPHERE_H

#include "Manifolds/Manifold.h"
#include "Manifolds/L2Sphere/L2SphereVariable.h"
#include "Manifolds/L2Sphere/L2SphereVector.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/* Riemannian Metric for the L2Sphere:
	TRAPEZOID: Composite trapezoidal rule*/
	enum Repa2NSMetric{ TRAPEZOID, L2SPHEREMETRICLENGTH };

	/*Retraction:
	L2SEXP: Exponential mapping
	NORMALIZED: R_x(etax) = (x + etax) / \|x + etax\|_{L^2}, this is not implemented here. */
	enum Repa2NSRetraction{ L2SEXP, NORMALIZED, L2SPHERERETRACTIONLENGTH };

	/*Vector transport:
	L2SEXP: parallel translation*/
	enum Repa2NSVectorTransport{ L2SPARALLELTRANSLATION, L2SPHEREVECTORTRANSPORTLENGTH };

	class L2Sphere : public Manifold{
	public:
		/*Construct the sphere in L^2([0, 1], R), inn uniformly-spaced points in [0, 1] are used to discretize functions.. */
		L2Sphere(integer inn);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~L2Sphere();

		/*The composite trapezoidal rule is used to define the Riemannian metric.*/
		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} v*/
		virtual void Projection(Variable *x, Vector *v, Vector *result) const;

		/*Exponential mapping is used*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*This is not done yet*/
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*This is not done yet*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Return 1; This manifold uses exponential mapping and parallel translation. Therefore, using beta = 1 satisfies the locking condition.*/
		virtual double Beta(Variable *x, Vector *etax) const;

		/*Parallel translation*/
		virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*Inverse parallel translation*/
		virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Compute result = H(:, start : end) * \mathcal{T}^{-1}, where H(:, start : end) denotes the matrix formed by columns from "start" to "end".
		\mathcal{T}^{-1} is the inverse parallel translation*/
		virtual void HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H(start : end, :), where H(start : end, :) denotes the matrix formed by cows from "start" to "end".
		\mathcal{T} is the parallel translation*/
		virtual void TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}. \mathcal{T} is the parallel translation*/
		virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

		/*Compute etaxflat = etax^{\flat}. */
		virtual void ObtainEtaxFlat(Variable *x, Vector *etax, Vector *etaxflat) const;

		/*This is not done yet*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*This is not done yet*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*result <-- v*/
		virtual void IntrProjection(Variable *x, Vector *v, Vector *result) const;

		/*etax is in the ambient space. This function projects etax onto the tangent space of x, i.e., result = P_{T_x M} v*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*When the metric is Euclidean, the Riemannian gradient is obtained by projecting the Euclidean
		onto the tangent space of x.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*The Riemannian action of the Hessian is obtained by Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;
	private:
		mutable integer n; /*The number of points to represent the continuous function*/

		Repa2NSMetric metric; /*Riemannian metric*/
		Repa2NSRetraction retraction; /*Retraction*/
		Repa2NSVectorTransport VecTran; /*Vector transport*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of L2SPHERE_H
