
#include "Manifolds/Grassmann/GrassVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	GrassVector::GrassVector(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	GrassVector *GrassVector::ConstructEmpty(void) const
	{
		return new GrassVector(size[0], size[1], size[2]);
	};
}; /*end of ROPTLIB namespace*/
