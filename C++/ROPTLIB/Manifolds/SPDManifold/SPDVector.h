/*
This file defines the class of a point on the tangent space of the manifold of symmetric positive definite matrices (SPD).

SmartSpace --> Element --> SPDVector

---- WH
*/

#ifndef SPDVECTOR_H
#define SPDVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDVector : public Element{
	public:
		/*Construct an empty vector on the tangent space of SPD. */
		SPDVector(integer row, integer col);

		/*Create an object of SPDVector with same size as this SPDVector.*/
		virtual SPDVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVECTOR_H
