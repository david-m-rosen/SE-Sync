
#include "Manifolds/CpxNStQOrth/CSOVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	CSOVariable::CSOVariable(integer r, integer c)
	{
		Element::Initialization(2, r * 2, c);
	};

	CSOVariable *CSOVariable::ConstructEmpty(void) const
	{
		return new CSOVariable(size[0] / 2, size[1]);
	};

	void CSOVariable::RandInManifold(void)
	{
		Element::RandGaussian();
	};
}; /*end of ROPTLIB namespace*/
