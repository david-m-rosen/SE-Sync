/*
 This file defines the class of a point on the preshape space.
 
 SmartSpace --> Element --> PSCVariable
 
 ---- YY
 */
#ifndef PSCVARIABLE_H
#define PSCVARIABLE_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"
#include "Problems/PreShapePathStraighten/PreShapePathStraighten.h"

/*Define the namespace*/
namespace ROPTLIB{
	class PSCVariable : public Element{
	public:
		// r is the number of points on a curve
		// l is the dimension of a space the curves are in
		// n is the number of curves on a path
    
		/*Construct an empty variable on the preshape space with only size information*/
		PSCVariable(integer r, integer l, integer n);
    
		/*Generate an initial path in the preshape space given the two end points*/
		void Generate(double *initial, double *end);
    
		/*Creat an object of PSCVariable with same size as this PSCVariable*/
		virtual PSCVariable *ConstructEmpty(void) const;
    
		/*This function randomly generates a point on the preshape manifold*/
		virtual void RandInManifold();
	};
}; /*end of ROPTLIB namespace*/

#endif // end of EUCVARIABLE_H
