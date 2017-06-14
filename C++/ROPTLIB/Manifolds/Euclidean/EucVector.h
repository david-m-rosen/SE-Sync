/*
This file defines the class of a point on the tangent space of the Euclidean space.

SmartSpace --> Element --> EucVector

---- WH
*/

#ifndef EUCVECTOR_H
#define EUCVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucVector : public Element{
	public:
		/*Construct an empty vector on the Euclidean space with only size information. */
		EucVector(integer r, integer c = 1, integer n = 1);

		/*Create an object of EucVector with same size as this EucVector.*/
		virtual EucVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVECTOR_H
