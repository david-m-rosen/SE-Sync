/*
This file defines the class of a point on the tangent space of the Stiefel manifold \St(p, n) = \{X \in R^{n \times p} | X^T X = I_p\}

SmartSpace --> Element --> StieVector

---- WH
*/

#ifndef STIEVECTOR_H
#define STIEVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class StieVector : public Element{
	public:
		/*Construct an empty vector on the tangent space of Stiefel manifold St(p, n) with only size information. */
		StieVector(integer n, integer p = 1, integer num = 1);

		/*Create an object of StieVector with same size as this StieVector.*/
		virtual StieVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of STIEVECTOR_H
