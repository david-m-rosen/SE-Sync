/*
This file defines the class of a point on the Oblique manifold Ob(n, num) = \{X \in R^{n times num} | diag(X^T X) = I_{num} \}.

SmartSpace --> ProductElement --> ObliqueVariable

---- WH
*/

#ifndef OBLIQUEVARIABLE_H
#define OBLIQUEVARIABLE_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/Sphere/SphereVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ObliqueVariable : public ProductElement{
	public:
		/*Construct an empty variable on the Oblique manifold with only size information.*/
		ObliqueVariable(integer n, integer num);

		/*Destruct by deleting all variables*/
		virtual ~ObliqueVariable();

		/*Create an object of ObliqueVariable with same size as this ObliqueVariable.*/
		virtual ObliqueVariable *ConstructEmpty(void) const;

		/*Print this point on oblique as a matrix*/
		virtual void Print(const char *name = "", bool isonlymain = true) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of OBLIQUEVARIABLE_H
