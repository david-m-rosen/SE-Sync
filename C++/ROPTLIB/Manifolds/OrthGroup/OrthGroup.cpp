
#include "Manifolds/OrthGroup/OrthGroup.h"

/*Define the namespace*/
namespace ROPTLIB{

	OrthGroup::OrthGroup(integer inn) :Stiefel(inn, inn)
	{
		name.assign("OrthGroup");
		delete EMPTYEXTR;
		delete EMPTYINTR;
		EMPTYEXTR = new OrthGroupVector(n, n);
		EMPTYINTR = new OrthGroupVector(IntrinsicDim);
	};

	OrthGroup::~OrthGroup(void)
	{
	};
}; /*end of ROPTLIB namespace*/
