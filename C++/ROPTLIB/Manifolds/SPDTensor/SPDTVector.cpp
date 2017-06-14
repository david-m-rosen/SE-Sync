
#include "Manifolds/SPDTensor/SPDTVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDTVector::SPDTVector(integer row, integer col, integer num)
	{
		SPDVector SPDV(row, col);
		Element **SPDVs = new Element *[num];
		for (integer i = 0; i < num; i++)
		{
			SPDVs[i] = &SPDV;
		}
		integer *powsintev = new integer[2];
		powsintev[0] = 0;
		powsintev[1] = num;

		ProductElementInitialization(SPDVs, num, powsintev, 1);

		delete[] powsintev;
		delete[] SPDVs;
	};

	SPDTVector *SPDTVector::ConstructEmpty(void) const
	{
		return new SPDTVector(elements[0]->Getsize()[0], elements[0]->Getsize()[1], numofelements);
	};
}; /*end of ROPTLIB namespace*/
