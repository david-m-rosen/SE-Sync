
#include "Manifolds/OrthGroup/OrthGroupVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	OrthGroupVariable::OrthGroupVariable(integer n) :StieVariable(n, n)
	{
	};

	OrthGroupVariable *OrthGroupVariable::ConstructEmpty(void) const
	{
		return new OrthGroupVariable(size[0]);
	};
}; /*end of ROPTLIB namespace*/
