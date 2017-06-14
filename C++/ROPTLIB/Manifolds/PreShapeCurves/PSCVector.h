/*
 This file defines the class of a point on the preshape space.
 
 SmartSpace --> Element --> PSCVector
 
 ---- YY
 */
#ifndef PSCVECTOR_H
#define PSCVECTOR_H

#include "Manifolds/Element.h"
#include <new>
#include <iostream>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{
	class PSCVector : public Element{
	public:
		/*Construct an empty vector on the preshape space with only the size information*/
		PSCVector(integer r, integer c, integer n);
    
		/*Creat an object of PSCVector with same size as this PSCVector*/
		virtual PSCVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCVECTOR_H
