/*
This file defines the class of a point on the manifold C_*^{n \times p} / O_p, where C_*^{n \times p} is a n by p full column rank
complex matrix and O_p is a p-by-p unitary group.

SmartSpace --> Element --> CSOVariable

---- WH
*/

#ifndef CSOVARIABLE_H
#define CSOVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

class CSOVariable : public Element{
	public:
		/*Construct an empty variable on the manifold C_*^{n \times p} / O_p with only size information. */
		CSOVariable(integer n, integer p = 1);

		/*Create an object of CSOVariable with same size as this CSOVariable.*/
		virtual CSOVariable *ConstructEmpty(void) const;

		/*This function randomly generates a point on the manifold.*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/

#endif // end of CSOVARIABLE_H
