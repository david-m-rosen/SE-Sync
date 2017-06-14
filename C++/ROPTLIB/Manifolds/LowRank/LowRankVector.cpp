
#include "Manifolds/LowRank/LowRankVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	LowRankVector::LowRankVector(integer Ur, integer Uc, integer Drc, integer Vr, integer Vc)
	{
		GrassVector dU(Ur, Uc);
		EucVector dD(Drc, Drc);
		GrassVector dV(Vr, Vc);

		Element **Elems = new Element *[3];
		Elems[0] = &dU;
		Elems[1] = &dD;
		Elems[2] = &dV;
		integer *powsintev = new integer[4];
		powsintev[0] = 0;
		powsintev[1] = 1;
		powsintev[2] = 2;
		powsintev[3] = 3;

		ProductElementInitialization(Elems, 3, powsintev, 3);

		delete[] powsintev;
		delete[] Elems;
	};

	LowRankVector::~LowRankVector(void)
	{
		/*If the sparse matrix exist in the temporary data, then invoke
		BLAS_usds to delete the sparse matrix. */
		if (TempDataExist("SparseMatrix"))
		{
			const SharedSpace *SM = this->ObtainReadTempData("SparseMatrix");
			/*If this sparse matrix is only referred by one element, then delete it.*/
			if (SM->GetSharedTimes()[0] == 1)
			{
				BLAS_usds(static_cast<blas_sparse_matrix> (SM->ObtainReadData()[0]));
			}
		}
	};

	LowRankVector *LowRankVector::ConstructEmpty(void) const
	{
		integer Ur = elements[0]->Getsize()[0];
		integer Uc = elements[0]->Getsize()[1];
		integer Drc = elements[1]->Getsize()[0];
		integer Vr = elements[2]->Getsize()[0];
		integer Vc = elements[2]->Getsize()[1];
		return new LowRankVector(Ur, Uc, Drc, Vr, Vc);
	};
}; /*end of ROPTLIB namespace*/
