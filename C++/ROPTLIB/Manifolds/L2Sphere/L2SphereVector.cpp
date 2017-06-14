
#include "Manifolds/L2Sphere/L2SphereVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	L2SphereVector::L2SphereVector(integer n)
	{
		Element::Initialization(1, n);
	};

	L2SphereVector *L2SphereVector::ConstructEmpty(void) const
	{
		return new L2SphereVector(size[0]);
	};
}; /*end of ROPTLIB namespace*/
