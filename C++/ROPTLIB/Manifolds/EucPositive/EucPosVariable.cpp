
#include "Manifolds/EucPositive/EucPosVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucPosVariable::EucPosVariable(integer r, integer c, integer n)
	{
		Element::Initialization(3, r, c, n);
	};

	EucPosVariable *EucPosVariable::ConstructEmpty(void) const
	{
		return new EucPosVariable(size[0], size[1], size[2]);
	};

	void EucPosVariable::RandInManifold(void)
	{
		Element::NewMemoryOnWrite();
		for (integer i = 0; i < length; i++)
			Space[i] = genrandreal();
	};
}; /*end of ROPTLIB namespace*/
