/*
This file defines the class of a point on the orthogonal group O_n = \{X \in R^{n \times n} | X^T X = I_n\}

SmartSpace --> Element --> StieVariable --> OrthGroupVariable

---- WH
*/

#ifndef ORTHGROUPVARIABLE_H
#define ORTHGROUPVARIABLE_H

#include "Manifolds/Stiefel/StieVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	class OrthGroupVariable : public StieVariable{
	public:
		/*Construct an empty variable on the orthogonal group O_n  with only size information. */
		OrthGroupVariable(integer n);

		/*Create an object of OrthGroupVariable with same size as this OrthGroupVariable.*/
		virtual OrthGroupVariable *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of ORTHGROUPVARIABLE_H