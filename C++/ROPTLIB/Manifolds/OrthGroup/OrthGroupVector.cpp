
#include "Manifolds/OrthGroup/OrthGroupVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	OrthGroupVector::OrthGroupVector(integer n, integer m) :StieVector(n, m)
	{
	};

	OrthGroupVector *OrthGroupVector::ConstructEmpty(void) const
	{
		return new OrthGroupVector(size[0], size[1]);
	};
}; /*end of ROPTLIB namespace*/
