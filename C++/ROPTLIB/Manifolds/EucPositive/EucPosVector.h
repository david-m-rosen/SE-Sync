/*
This file defines the class of a point on the tangent cone of the EucPositive space.

SmartSpace --> Element --> EucPosVector

---- WH
*/

#ifndef EUCPOSVECTOR_H
#define EUCPOSVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucPosVector : public Element{
	public:
		/*Construct an empty vector on the EucPostive space with only size information. */
		EucPosVector(integer r, integer c = 1, integer n = 1);

		/*Create an object of EucPosVector with same size as this EucPosVector.*/
		virtual EucPosVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCPOSVECTOR_H
