
#include "Manifolds/SPDManifold/SPDVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDVariable::SPDVariable(integer n)
	{
		Element::Initialization(2, n, n);
	};

	SPDVariable *SPDVariable::ConstructEmpty(void) const
	{
		return new SPDVariable(size[0]);
	};

	void SPDVariable::RandInManifold(void)
	{
		integer n = size[0];
		/*temp is an n by n matrix*/
		double *temp = new double[n * n];
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i; j < n; j++)
			{
				temp[i + j * n] = 0;
				temp[j + i * n] = genrandnormal();
			}
		}

		NewMemoryOnWrite();
		/*Space <-- temp * temp^T. Thus, Space points to a symmetric positive definite matrix. Therefore,
		this SPDVariable is a SPD matrix.*/
		dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, temp, &n, temp, &n, &GLOBAL::DZERO, Space, &n);
		delete[] temp;
	};
}; /*end of ROPTLIB namespace*/
