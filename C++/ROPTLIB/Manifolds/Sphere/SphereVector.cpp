
#include "Manifolds/Sphere/SphereVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereVector::SphereVector(integer n) :StieVector(n)
	{
	};

	SphereVector *SphereVector::ConstructEmpty(void) const
	{
		return new SphereVector(length);
	};
}; /*end of ROPTLIB namespace*/
