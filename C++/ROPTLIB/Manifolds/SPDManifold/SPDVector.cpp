
#include "Manifolds/SPDManifold/SPDVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDVector::SPDVector(integer row, integer col)
	{
		Element::Initialization(2, row, col);
	};

	SPDVector *SPDVector::ConstructEmpty(void) const
	{
		return new SPDVector(size[0], size[1]);
	};
}; /*end of ROPTLIB namespace*/
