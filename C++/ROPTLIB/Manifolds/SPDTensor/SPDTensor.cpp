
#include "Manifolds/SPDTensor/SPDTensor.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDTensor::SPDTensor(integer indim, integer innum) : ProductManifold(1, new SPDManifold(indim), innum)
	{
		name.assign("SPDTensor");
		dim = indim;
		num = innum;
		delete EMPTYEXTR;
		delete EMPTYINTR;
		EMPTYEXTR = new SPDTVector(indim, indim, innum);
		EMPTYINTR = new SPDTVector(indim * (indim + 1) / 2, 1, innum);
	};

	SPDTensor::~SPDTensor(void)
	{
		for (integer i = 0; i < numofmani; i++)
		{
			delete manifolds[i];
		}
	};

	void SPDTensor::CholeskyRepresentation(Variable *x) const
	{
		const SharedSpace * Xalpha = x->ObtainReadTempData("Xalpha");
		const double *Xalphaptr = Xalpha->ObtainReadData();
		integer N = Xalpha->Getsize()[2];

		SharedSpace *SharedL = new SharedSpace(3, dim, dim, N);
		double *LM = SharedL->ObtainWriteEntireData();
		for (integer k = 0; k < N; k++)
		{
			for (integer i = 0; i < dim; i++)
			{
				for (integer j = i; j < dim; j++)
				{
					LM[i + j * dim + k * dim * dim] = 0;
					LM[j + i * dim + k * dim * dim] = Xalphaptr[j + i * dim + k * dim * dim];
				}
			}
		}

		integer info, ddim = dim;
		for (integer k = 0; k < N; k++)
		{
			dpotrf_(GLOBAL::L, &ddim, LM + k * dim * dim, &ddim, &info);
			if (info != 0)
			{
				printf("Warning: SPDTensor::CholeskyRepresentation fails with info:%d!\n", info);
			}
		}
		x->AddToTempData("XaL", SharedL);
	};
}; /*end of ROPTLIB namespace*/
