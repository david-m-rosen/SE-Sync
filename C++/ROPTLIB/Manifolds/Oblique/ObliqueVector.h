/*
This file defines the class of a point on the tangent space of the Oblique manifold 
Ob(n, num) = \{X \in R^{n times num} | diag(X^T X) = I_{num} \}.

SmartSpace --> ProductElement --> ObliqueVector

---- WH
*/

#ifndef OBLIQUEVECTOR_H
#define OBLIQUEVECTOR_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/Sphere/SphereVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ObliqueVector : public ProductElement{
	public:
		/*Construct an empty vector on the tangent space of Oblique manifold with only size information. */
		ObliqueVector(integer n, integer num);

		/*Destruct by deleting all components*/
		virtual ~ObliqueVector();

		/*Create an object of ObliqueVector with same size as this ObliqueVector.*/
		virtual ObliqueVector *ConstructEmpty(void) const;

		/*Print this point on tangent space of oblique as a matrix*/
		virtual void Print(const char *name = "", bool isonlymain = true) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of OBLIQUEVECTOR_H
