/*
This file defines the class of a point on the Stiefel manifold \St(p, n) = \{X \in R^{n \times p} | X^T X = I_p\}

SmartSpace --> Element --> StieVariable

---- WH
*/

#ifndef STIEVARIABLE_H
#define STIEVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class StieVariable : public Element{
	public:
		/*Construct an empty variable on the Stiefel manifold St(p, n) with only size information. */
		StieVariable(integer n, integer p = 1, integer num = 1);

		/*Create an object of StieVariable with same size as this StieVariable.*/
		virtual StieVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the Stiefel manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/

#endif // end of STIEVARIABLE_H
