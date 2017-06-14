
#include "Manifolds/Euclidean/EucVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucVariable::EucVariable(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	EucVariable *EucVariable::ConstructEmpty(void) const
	{
		return new EucVariable(size[0], size[1], size[2]);
	};

	void EucVariable::RandInManifold(void)
	{
		Element::RandGaussian();
	};
}; /*end of ROPTLIB namespace*/
