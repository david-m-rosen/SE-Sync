
#include "Manifolds/Oblique/ObliqueVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	ObliqueVector::ObliqueVector(integer n, integer num)
	{
		SphereVector SV(n);

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

	ObliqueVector::~ObliqueVector(void)
	{
	};

	ObliqueVector *ObliqueVector::ConstructEmpty(void) const
	{
		return new ObliqueVector(elements[0]->Getlength(), numofelements);
	};

	void ObliqueVector::Print(const char *name, bool isonlymain) const
	{
		if (isonlymain)
		{
			if (Space == nullptr)
			{
				if (size == nullptr)
				{
					printf("%s is an empty data with size 0", name);
				}
				else
				{
					printf("%s is an empty data with size %d", name, size[0]);
				}
				for (integer i = 1; i < ls; i++)
					printf(" x %d", size[i]);
				printf("\n");
				return;
			}
			printf("%s , shared times:%d, shared times address:%p\n", name, *sharedtimes, sharedtimes);
			integer n = elements[0]->Getlength();
			integer num = numofelements;
			for (integer i = 0; i < n; i++)
			{
				for (integer j = 0; j < num; j++)
				{
					printf("%f\t", elements[j]->GetSpace()[i]);
				}
				printf("\n");
			}
			return;
		}
		ProductElement::Print(name, isonlymain);
	};
}; /*end of ROPTLIB namespace*/
