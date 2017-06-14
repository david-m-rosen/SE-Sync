
#include "Manifolds/CpxNStQOrth/CSOVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	CSOVector::CSOVector(integer r, integer c)
	{
		Element::Initialization(2, r, c);
	};

	CSOVector *CSOVector::ConstructEmpty(void) const
	{
		return new CSOVector(size[0], size[1]);
	};
}; /*end of ROPTLIB namespace*/
