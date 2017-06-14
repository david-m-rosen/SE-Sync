
#include "Manifolds/Sphere/SphereVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereVariable::SphereVariable(integer n) :StieVariable(n)
	{
	};

	SphereVariable *SphereVariable::ConstructEmpty(void) const
	{
		return new SphereVariable(length);
	};
}; /*end of ROPTLIB namespace*/
