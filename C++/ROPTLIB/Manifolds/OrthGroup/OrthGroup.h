/*
This file defines the class for the orthogonal group O_n = \{X \in R^{n \times n} | X^T X = I_n\}
It defines the common properties and features of the manifold.

Manifold --> Stiefel --> OrthGroup

---- WH
*/

#ifndef ORTHGROUP_H
#define ORTHGROUP_H

#include "Manifolds/OrthGroup/OrthGroupVariable.h"
#include "Manifolds/OrthGroup/OrthGroupVector.h"
#include "Manifolds/Stiefel/Stiefel.h"

/*Define the namespace*/
namespace ROPTLIB{

	class OrthGroup : public Stiefel{
	public:
		/*Construct the orthogonal group with size inn by inn*/
		OrthGroup(integer inn);

		/*Delete the orthogonal group*/
		virtual ~OrthGroup();
	};
}; /*end of ROPTLIB namespace*/
#endif // end of ORTHGROUP_H
