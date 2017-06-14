/*
This file defines the class of a point on the product of symmetric positive definite matrices manifold (SPDTensor).

SmartSpace --> Element --> SPDTVariable

---- WH
*/

#ifndef SPDTVARIABLE_H
#define SPDTVARIABLE_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/SPDManifold/SPDVariable.h"
#include <new>
#include <iostream>
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDTVariable : public ProductElement{
	public:
		/*Construct an empty variable on SPDTensor with only size information. */
		SPDTVariable(integer dim, integer num);

		/*Create an object of SPDTVariable with same size as this SPDTVariable.*/
		virtual SPDTVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCTVARIABLE_H
