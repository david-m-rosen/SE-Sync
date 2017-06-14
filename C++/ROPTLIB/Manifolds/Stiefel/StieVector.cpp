
#include "Manifolds/Stiefel/StieVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	StieVector::StieVector(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	StieVector *StieVector::ConstructEmpty(void) const
	{
		return new StieVector(size[0], size[1], size[2]);
	};
}; /*end of ROPTLIB namespace*/
