
#include "Manifolds/Oblique/Oblique.h"

/*Define the namespace*/
namespace ROPTLIB{

	Oblique::Oblique(integer n, integer num) : ProductManifold(1, new Sphere(n), num)
	{
		name.assign("Oblique");
		delete EMPTYEXTR;
		delete EMPTYINTR;
		EMPTYEXTR = new ObliqueVector(n, num);
		EMPTYINTR = new ObliqueVector(n - 1, num);
	};

	Oblique::~Oblique(void)
	{
		for (integer i = 0; i < numofmani; i++)
		{
			delete manifolds[i];
		}
	};

	void Oblique::ChooseObliqueParamsSet1(void)
	{
		Sphere *S = dynamic_cast<Sphere *> (manifolds[0]);
		S->ChooseStieParamsSet1();
	};

	void Oblique::ChooseObliqueParamsSet2(void)
	{
		Sphere *S = dynamic_cast<Sphere *> (manifolds[0]);
		S->ChooseSphereParamsSet2();
		integer n = S->GetExtrDim();
		integer num = numoftotalmani;
		delete EMPTYINTR;
		EMPTYINTR = new ObliqueVector(n, num);
		//SetEMPTYINTR();
	};

	void Oblique::ChooseObliqueParamsSet3(void)
	{
		Sphere *S = dynamic_cast<Sphere *> (manifolds[0]);
		S->ChooseSphereParamsSet3();
		integer n = S->GetExtrDim();
		integer num = numoftotalmani;
		delete EMPTYINTR;
		EMPTYINTR = new ObliqueVector(n, num);
		//SetEMPTYINTR();
	};

	void Oblique::ChooseObliqueParamsSet4(void)
	{
		Sphere *S = dynamic_cast<Sphere *> (manifolds[0]);
		S->ChooseSphereParamsSet4();
		integer n = S->GetExtrDim();
		integer num = numoftotalmani;
		delete EMPTYINTR;
		EMPTYINTR = new ObliqueVector(n, num);
		//SetEMPTYINTR();
	};

	void Oblique::ChooseObliqueParamsSet5(void)
	{
		Sphere *S = dynamic_cast<Sphere *> (manifolds[0]);
		S->ChooseSphereParamsSet5();
		//SetEMPTYINTR();
	};

	void Oblique::SetParams(PARAMSMAP params)
	{
		Manifold::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("ParamSet"))
			{
				switch (static_cast<integer> (iter->second))
				{
				case 1:
					ChooseObliqueParamsSet1();
					break;
				case 2:
					ChooseObliqueParamsSet2();
					break;
				case 3:
					ChooseObliqueParamsSet3();
					break;
				case 4:
					ChooseObliqueParamsSet4();
					break;
				case 5:
					ChooseObliqueParamsSet5();
					break;
				default:
					break;
				}
			}
		}
	};
}; /*end of ROPTLIB namespace*/
