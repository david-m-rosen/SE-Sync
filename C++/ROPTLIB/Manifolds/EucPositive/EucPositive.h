/*
This file defines the class for the Eucldean space with all entries positive.
Note that this is not a manifold. This class defines objects for gradient projection
methods.

Manifold --> EucPositive

---- WH
*/
#ifndef EUCPOSITIVE_H
#define EUCPOSITIVE_H


#include "Manifolds/EucPositive/EucPosVariable.h"
#include "Manifolds/EucPositive/EucPosVector.h"
#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucPositive : public Manifold{
	public:
		/*Construct the EucPositive space*/
		EucPositive(integer r, integer c = 1, integer n = 1);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~EucPositive(void);

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*Project etax onto the tangent cone at x.*/
		virtual void Projection(Variable *x, Vector *etax, Vector *result) const;

		/*Compute the retraction result = R_x(etax)
		Default: result = P (x + etax);*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*We project the Euclidean gradient to the tangent cone at $x$, i.e., gf <-- P egf*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*TODO: xix <-- exix*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		integer row; /*The first dimension of the space, i.e., the number of rows */
		integer col; /*The second dimension of the space, i.e., the number of columns */
		integer num; /*The third dimension of the space*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCPOSITIVE_H
