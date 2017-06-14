
#include "Manifolds/PreShapeCurves/PSCVector.h"

/*Define the namespace*/
namespace ROPTLIB{
	PSCVector::PSCVector(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};
    
	PSCVector *PSCVector::ConstructEmpty(void) const
	{
		return new PSCVector(size[0], size[1], size[2]);
	};
}; /*end of ROPTLIB namespace*/
