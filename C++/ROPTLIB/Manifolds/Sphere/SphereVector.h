/*
This file defines the class of a point on the tangent space of the unit Sphere S^{n-1} = \{x \in R^{n} | x^T x = 1\}

SmartSpace --> Element --> StieVector --> SphereVector

---- WH
*/

#ifndef SPHEREVECTOR_H
#define SPHEREVECTOR_H

#include "Manifolds/Stiefel/StieVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereVector :public StieVector{
	public:
		/*Construct an empty vector on the tangent space of the unit sphere S^{n-1} with only size information. */
		SphereVector(integer n);

		/*Create an object of SphereVector with same size as this SphereVector.*/
		virtual SphereVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of SPHEREVECTOR_H