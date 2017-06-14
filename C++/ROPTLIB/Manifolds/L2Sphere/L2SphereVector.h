/*
This file defines the class of a point on the tangent space of the sphere in L^2([0, 1], R).

SmartSpace --> Element --> EucVector

---- WH
*/

#ifndef L2SPHEREVECTOR_H
#define L2SPHEREVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class L2SphereVector : public Element{
	public:
		/*Construct an empty vector on the sphere of L^2([0, 1], R) with only size information.
		n denotes the number of points to represent the continuous function.*/
		L2SphereVector(integer n);

		/*Create an object of L2SphereVector with same size as this L2SphereVector.*/
		virtual L2SphereVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of L2SPHEREVECTOR_H
