/*
This file defines the class of a point on the tangent space of the Grassmann manifold \Gr(p, n) = \{[X]| X^T X = I_p, X \in R^{n \times p}\}
and [X] = \{XO | O^T O = I_p, O \in R^{p \times p}\}

SmartSpace --> Element --> GrassVector

---- WH
*/

#ifndef GRASSVECTOR_H
#define GRASSVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class GrassVector : public Element{
	public:
		/*Construct an empty vector on the tangent space of Grassmann manifold Gr(p, n) with only size information. */
		GrassVector(integer n, integer p = 1, integer num = 1);

		/*Create an object of GrassVector with same size as this GrassVector.*/
		virtual GrassVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of GRASSVECTOR_H
