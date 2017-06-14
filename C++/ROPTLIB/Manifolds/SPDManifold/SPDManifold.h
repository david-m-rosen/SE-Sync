/*
This file defines the class for the manifold of symmetric positive definite matrices (SPDManifold).
The affine invariant metric is used. Intrinsic representation of tangent vectors are used.
See details in: 
	[YHAG15] Xinru Yuan, Wen Huang, P.-A. Absil, K. A. Gallivan. "A Riemannian Limited-memory BFGS 
	Algorithm for Computing the Matrix Geometric Mean".

Manifold --> SPD

---- WH
*/
#ifndef SPDMANIFOLD_H
#define SPDMANIFOLD_H

#include "Manifolds/SPDManifold/SPDVariable.h"
#include "Manifolds/SPDManifold/SPDVector.h"
#include "Manifolds/Manifold.h"
#include "Others/MyMatrix.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDManifold : public Manifold{
	public:
		/*Construct the SPD manifold*/
		SPDManifold(integer inn);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~SPDManifold(void);

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*The second order approximation of the exponential mapping, i.e., result = x + etax + 0.5 etax x^{-1} etax */
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*computes beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		Default: beta <-- 1*/
		virtual double Beta(Variable *x, Vector *etax) const;

		/*compute the distance between two points on the manifold.
		Default: the distance under affine invariant metric: ||log(x1^{-1/2) x2 x1^{-1/2}||_F*/
		virtual double Dist(Variable *x1, Variable *x2) const;

		/*The vector transport by differentiated the retraction, which is R_x(etax) = x + etax + 0.5 etax x^{-1} etax.
		Therefore, it is xix + 0.5 xix x^{-1} etax + 0.5 etax x^{-1} xix */
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*TODO*/
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*The orthogonal projection onto the tangent space using the affine invariance metric is: result = (etax + etax^T) / 2*/
		virtual void ExtrProjection(Variable *x, Vector *etax, Vector *result) const;

		/*Compute the intrinsic representation of a tagnent vector etax, intreta = upper triangle
		of L^{-1} etax L^{-T}, where x = L L^T. */
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Compute the extrinsic representation of a tagnent vector intretax. Inverse operation of the function ObtainIntr */
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*gf = x * egf * x */
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*Not done yet. Temporarily use: xix <-- exix*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*X = L L^T, the L is attached to X as a temporary data.*/
		void CholeskyRepresentation(Variable *x) const;

	protected:
		integer n; /*The size of the space, i.e., the number of row/column */
	};
}; /*end of ROPTLIB namespace*/
#endif // end of SPDMANIFOLD_H
