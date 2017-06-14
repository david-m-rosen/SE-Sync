/*
This file defines the class of a point on the unit Sphere S^{n-1} = \{x \in R^{n} | x^T x = 1\}

SmartSpace --> Element --> StieVariable --> SphereVariable

---- WH
*/

#ifndef SPHEREVARIABLE_H
#define SPHEREVARIABLE_H

#include "Manifolds/Stiefel/StieVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SphereVariable : public StieVariable{
	public:
		/*Construct an empty variable on the unit Sphere S^{n-1} with only size information. */
		SphereVariable(integer n);

		/*Create an object of SphereVariable with same size as this SphereVariable.*/
		virtual SphereVariable *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of SPHEREVARIABLE_H