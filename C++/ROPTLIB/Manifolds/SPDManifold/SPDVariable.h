/*
This file defines the class of a point on the the manifold of symmetric positive definite matrices (SPD)

SmartSpace --> Element --> SPDVariable

---- WH
*/

#ifndef SPDVARIABLE_H
#define SPDVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDVariable : public Element{
	public:
		/*Construct an empty variable on SPD with only size information. */
		SPDVariable(integer n);

		/*Create an object of SPDVariable with same size as this SPDVariable.*/
		virtual SPDVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVARIABLE_H
