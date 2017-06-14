/*
This file defines the class of a point on the tangent space of the product of symmetric positive definite matrices manifold (SPDTensor).

SmartSpace --> Element --> SPDTVector

---- WH
*/

#ifndef SPDTVECTOR_H
#define SPDTVECTOR_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/SPDManifold/SPDVector.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDTVector : public ProductElement{
	public:
		/*Construct an empty vector on the tangent space of SPDTensor. */
		SPDTVector(integer row, integer col, integer num);

		/*Create an object of SPDTVector with same size as this SPDTVector.*/
		virtual SPDTVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVECTOR_H
