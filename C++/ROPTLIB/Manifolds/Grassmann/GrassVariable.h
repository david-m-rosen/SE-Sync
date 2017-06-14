/*
This file defines the class of a point on the Grassmann manifold \Gr(p, n) = \{[X]| X^T X = I_p, X \in R^{n \times p}\}
and [X] = \{XO | O^T O = I_p, O \in R^{p \times p}\}

SmartSpace --> Element --> GrassVariable

---- WH
*/

#ifndef GRASSVARIABLE_H
#define GRASSVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class GrassVariable : public Element{
	public:
		/*Construct an empty variable on the Grassmann manifold Gr(p, n) with only size information. */
		GrassVariable(integer n, integer p = 1, integer num = 1);

		/*Create an object of GrassVariable with same size as this GrassVariable.*/
		virtual GrassVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the Grassmann manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of GRASSVARIABLE_H
