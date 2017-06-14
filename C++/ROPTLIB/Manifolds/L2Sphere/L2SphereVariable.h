/*
This file defines the class of a point on the sphere in L^2([0, 1], R) represented by its values on uniformly-spaced
grid, i.e., ( x(0/(n-1)), x(1/(n-1)), cdots, x((n-1)/(n-1)) ).

SmartSpace --> Element --> L2SphereVariable

---- WH
*/

#ifndef L2SPHEREVARIABLE_H
#define L2SPHEREVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class L2SphereVariable : public Element{
	public:
		/*Construct an empty variable on the Euclidean space with only size information.
		n denotes the number of points to represent the continuous function.*/
		L2SphereVariable(integer n);

		/*Create an object of L2SphereVariable with same size as this L2SphereVariable.*/
		virtual L2SphereVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of L2SPHEREVARIABLE_H
