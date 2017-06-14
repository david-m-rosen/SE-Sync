/*
This file defines the class for storing a linear operator. 

SmartSpace --> LinearOPE

---- WH
*/

#ifndef LINEAROPE_H
#define LINEAROPE_H

#include "Manifolds/Element.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LinearOPE : public SmartSpace{
	public:
		/*Defines an linear operator on R^{s times s}. This is an empty operator with only size information*/
		LinearOPE(integer s);

		/*Create an object of LinearOPE with same size as this LinearOPE.*/
		virtual LinearOPE *ConstructEmpty(void) const;

		/*Assign this LinearOPE to be scalar * identity */
		virtual void ScaledIdOPE(double scalar = 1);
	};
}; /*end of ROPTLIB namespace*/
#endif // end of LINEAROPE_H
