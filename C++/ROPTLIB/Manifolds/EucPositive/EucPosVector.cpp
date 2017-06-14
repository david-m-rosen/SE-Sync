
#include "Manifolds/EucPositive/EucPosVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucPosVector::EucPosVector(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	EucPosVector *EucPosVector::ConstructEmpty(void) const
	{
		return new EucPosVector(size[0], size[1], size[2]);
	};
}; /*end of ROPTLIB namespace*/
