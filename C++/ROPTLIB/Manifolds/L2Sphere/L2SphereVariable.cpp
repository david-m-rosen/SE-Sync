
#include "Manifolds/L2Sphere/L2SphereVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	L2SphereVariable::L2SphereVariable(integer n)
	{
		Element::Initialization(1, n);
	};

	L2SphereVariable *L2SphereVariable::ConstructEmpty(void) const
	{
		return new L2SphereVariable(size[0]);
	};

	void L2SphereVariable::RandInManifold(void)
	{
		this->RandGaussian();
		double norm = Space[0] * Space[0] / 2;
		for (integer i = 1; i < length - 1; i++)
		{
			norm += Space[i] * Space[i];
		}
		norm += Space[length - 1] * Space[length - 1] / 2;
		norm /= (length - 1);
		norm = sqrt(norm);
		double a = 1.0 / norm;
		integer inc = 1;
		// Space <- a * Space, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
		dscal_(&length, &a, Space, &inc);
	};
}; /*end of ROPTLIB namespace*/
