/*
This file defines the class of a point on the Euclidean space

SmartSpace --> Element --> EucVariable

---- WH
*/

#ifndef EUCVARIABLE_H
#define EUCVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucVariable : public Element{
	public:
		/*Construct an empty variable on the Euclidean space with only size information. */
		EucVariable(integer r, integer l = 1, integer n = 1);

		/*Create an object of EucVariable with same size as this EucVariable.*/
		virtual EucVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVARIABLE_H
