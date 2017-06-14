
#include "Manifolds/Euclidean/EucVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucVector::EucVector(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	EucVector *EucVector::ConstructEmpty(void) const
	{
		return new EucVector(size[0], size[1], size[2]);
	};
}; /*end of ROPTLIB namespace*/
