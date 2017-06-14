
#include "Manifolds/SPDTensor/SPDTVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDTVariable::SPDTVariable(integer dim, integer num)
	{
		SPDVariable SV(dim);

		Element **SVs = new Element *[num];
		for (integer i = 0; i < num; i++)
		{
			SVs[i] = &SV;
		}
		integer *powsintev = new integer[2];
		powsintev[0] = 0;
		powsintev[1] = num;

		ProductElementInitialization(SVs, num, powsintev, 1);

		delete[] powsintev;
		delete[] SVs;
	};

	SPDTVariable *SPDTVariable::ConstructEmpty(void) const
	{
		return new SPDTVariable(elements[0]->Getsize()[0], numofelements);
	};

	void SPDTVariable::RandInManifold(void)
	{
		integer n = elements[0]->Getsize()[0], num = numofelements;
		/*temp is an n by n matrix*/
		double *temp = new double[n * n];
		NewMemoryOnWrite();

		for (integer k = 0; k < num; k++)
		{
			for (integer i = 0; i < n; i++)
			{
				for (integer j = 0; j < n; j++)
				{
					temp[j + i * n] = genrandnormal();
				}
			}
			/*Space <-- temp * temp^T. Thus, Space points to a symmetric positive definite matrix. Therefore,
			this SPDVariable is a SPD matrix.*/
			dgemm_(GLOBAL::N, GLOBAL::T, &n, &n, &n, &GLOBAL::DONE, temp, &n, temp, &n, &GLOBAL::DZERO, Space + k * n * n, &n);
		}

		delete[] temp;
	};
}; /*end of ROPTLIB namespace*/
