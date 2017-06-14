
#include "Solvers/RGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	RGS::RGS(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void RGS::SetDefaultParams(void)
	{
		SolversLSLPSub::SetDefaultParams();
		Xs = nullptr;
		InitSteptype = ONESTEP;
		LineSearch_LS = ARMIJO;
		LS_ratio1 = 0.25;
		LS_ratio2 = 0.25;
		SolverName.assign("RGS");
	};

	RGS::~RGS(void)
	{
		DeleteVectors(Xs, Lengthgfs);
	};

	void RGS::GetSearchDir(void)
	{
		gf1->CopyTo(gfs[0]);
		Mani->RandomTangentVectors(x1, Lengthgfs - 1, gfs + 1);
		
		/*normalized all the tangent vectors such that their norms equal Eps*/
		double tmp = 0;
		for (integer i = 1; i < Lengthgfs; i++)
		{
			tmp = sqrt(Mani->Metric(x1, gfs[i], gfs[i]));
			Mani->ScaleTimesVector(x1, Eps / tmp, gfs[i], gfs[i]);
		}

		/*Apply retraction and obtain points around x1*/
		/*Compute gradients at Xs[i] and transport them to the tangent space at x1*/
		for (integer i = 1; i < Lengthgfs; i++)
		{
			Mani->Retraction(x1, gfs[i], Xs[i]); nR++;
			Prob->f(Xs[i]); nf++;
			Prob->Grad(Xs[i], zeta); ng++;
			Mani->InverseVectorTransport(x1, gfs[i], Xs[i], zeta, zeta);
			zeta->CopyTo(gfs[i]);
		}

		ngf = sqrt(MinPNormConHull(Mani, x1, gfs, Lengthgfs, nullptr, nullptr, gf, nullptr, 0));
		subprobtimes++;

		/*eta1 is viewed as the search direction*/
		Mani->ScaleTimesVector(x1, -1.0 / ngf, gf, eta1);
	};
}; /*end of ROPTLIB namespace*/
