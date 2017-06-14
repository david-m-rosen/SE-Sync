
#include "Manifolds/LowRank/LowRankVariable.h"

/*Define the namespace*/
namespace ROPTLIB{

	LowRankVariable::LowRankVariable(integer m, integer n, integer r)
	{
		GrassVariable U(m, r);
		EucVariable D(r, r);
		GrassVariable V(n, r);

		Element **Elems = new Element *[3];
		Elems[0] = &U;
		Elems[1] = &D;
		Elems[2] = &V;
		integer *powsintev = new integer[4];
		powsintev[0] = 0;
		powsintev[1] = 1;
		powsintev[2] = 2;
		powsintev[3] = 3;

		ProductElementInitialization(Elems, 3, powsintev, 3);

		delete[] powsintev;
		delete[] Elems;
	};

	LowRankVariable::~LowRankVariable(void)
	{
	};

	LowRankVariable *LowRankVariable::ConstructEmpty(void) const
	{
		integer m = elements[0]->Getsize()[0];
		integer r = elements[1]->Getsize()[0];
		integer n = elements[2]->Getsize()[0];
		return new LowRankVariable(m, n, r);
	};
}; /*end of ROPTLIB namespace*/
