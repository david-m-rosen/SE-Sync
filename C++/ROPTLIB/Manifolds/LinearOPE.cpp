
#include "Manifolds/LinearOPE.h"

/*Define the namespace*/
namespace ROPTLIB{

	LinearOPE::LinearOPE(integer s)
	{
		SmartSpace::Initialization(2, s, s);
	};

	LinearOPE *LinearOPE::ConstructEmpty(void) const
	{
		return new LinearOPE(*size);
	};

	void LinearOPE::ScaledIdOPE(double scalar)
	{
		SmartSpace::NewMemoryOnWrite();
		integer ell = size[0];
		for (integer i = 0; i < ell; i++)
		{
			Space[i + i * ell] = scalar;
			for (integer j = i + 1; j < ell; j++)
			{
				Space[i + j * ell] = 0;
				Space[j + i * ell] = 0;
			}
		}
	};
}; /*end of ROPTLIB namespace*/
