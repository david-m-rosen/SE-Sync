/*
This file defines the class of a point on the EucPositive space

SmartSpace --> Element --> EucPosVariable

---- WH
*/

#ifndef EUCPOSVARIABLE_H
#define EUCPOSVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucPosVariable : public Element{
	public:
		/*Construct an empty variable on the EucPositive space with only size information. */
		EucPosVariable(integer r, integer l = 1, integer n = 1);

		/*Create an object of EucPosVariable with same size as this EucPosVariable.*/
		virtual EucPosVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVARIABLE_H
