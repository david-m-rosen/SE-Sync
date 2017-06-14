/*
This file defines the class of a point on the tangent space of the orthogonal group O_n = \{X \in R^{n \times n} | X^T X = I_n\}

SmartSpace --> Element --> StieVector --> OrthGroupVector

---- WH
*/

#ifndef ORTHGROUPVECTOR_H
#define ORTHGROUPVECTOR_H

#include "Manifolds/Stiefel/StieVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	class OrthGroupVector :public StieVector{
	public:
		/*Construct an empty vector on the tangent space of orthogonal group O_n with only size information. */
		OrthGroupVector(integer n, integer m = 1);

		/*Create an object of OrthGroupVector with same size as this OrthGroupVector.*/
		virtual OrthGroupVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of ORTHGROUPVECTOR_H
