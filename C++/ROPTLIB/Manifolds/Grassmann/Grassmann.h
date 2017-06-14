/*
This file defines the class for the Grassmann manifold \Gr(p, n) = \{[X]| X^T X = I_p, X \in R^{n \times p}\}
and [X] = \{XO | O^T O = I_p, O \in R^{p \times p}\}.
It defines the common properties and features of the manifold.

Manifold --> Grassmann

---- WH
*/

#ifndef GRASSMANN_H
#define GRASSMANN_H

//#include <cmath>

#include "Manifolds/Grassmann/GrassVariable.h"
#include "Manifolds/Grassmann/GrassVector.h"
#include "Manifolds/Manifold.h"
#include "Others/MyMatrix.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Grassmann : public Manifold{
	public:
		/*Construct the Grassmann manifold: Gr(p, n) and set up default parameters*/
		Grassmann(integer n, integer p);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~Grassmann(void);

		/*Compute the qf retraction defined in [AMS2008, (4.8)].
			[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
			Princeton University Press, Princeton, NJ, 2008.*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*This cotangent vector is used in the RBFGS defined in [RW2012]
			[RW2012]: W. Ring and B. Wirth. Optimization methods on Riemannian manifolds and their application to shape space.
			SIAM Journal on Optimization, 22(2):596?27, January 2012. */
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*The vector transport by differentiated the retraction of the qf retraction */
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Obtain beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		beta has computed in "DiffRetraction". It is not necessary to recompute it in this function. */
		virtual double Beta(Variable *x, Vector *etax) const;

		/*Call a member function "ObtainIntrHHR" or "ObtainIntrSquare" based on member variable "retraction". */
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Call a member function "ObtainExtrHHR" or "ObtainExtrSquare" based on member variable "retraction". */
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*ExtrProjection is given by: result = v - x x^T v.*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*the Riemannian gradient is obtained by projecting the Euclidean onto the tangent space of x.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*From the Euclidean action of the Hessian exix:=D grad f[etax] to the Riemannian action of the Hessian:
		xix:=\nabla_etax Rgrad f(x) = (I - x x^T) (exix - etax x^T Egradf)*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

	protected:
		integer n; /*The number of row*/
		integer p; /*The number of column*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of GRASSMANN_H
