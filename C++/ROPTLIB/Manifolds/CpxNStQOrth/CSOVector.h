/*
This file defines the class of a point on the tangent space of C_*^{n \times p} / O_p, where C_*^{n \times p} is a n by p full column rank
complex matrix and O_p is a p-by-p unitary group.

SmartSpace --> Element --> CSOVector

---- WH
*/

#ifndef CSOVECTOR_H
#define CSOVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class CSOVector : public Element{
	public:
		/*Construct an empty variable on the tangent space of C_*^{n \times p} / O_p with only size information. */
		CSOVector(integer n, integer p = 1);

		/*Create an object of CSOVector with same size as this CSOVector.*/
		virtual CSOVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of EUCVECTOR_H
