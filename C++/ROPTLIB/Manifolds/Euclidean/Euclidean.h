/*
This file defines the class for the Eucldean space.

Manifold --> Euclidean

---- WH
*/
#ifndef EUCLIDEAN_H
#define EUCLIDEAN_H


#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Euclidean : public Manifold{
	public:
		/*Construct the Euclidean space*/
		Euclidean(integer r, integer c = 1, integer n = 1);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~Euclidean(void);

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*gf <-- egf*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*xix <-- exix*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		integer row; /*The first dimension of the space, i.e., the number of rows */
		integer col; /*The second dimension of the space, i.e., the number of columns */
		integer num; /*The third dimension of the space*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCLIDEAN_H
