
#include "Solvers/RBroydenFamily.h"

/*Define the namespace*/
namespace ROPTLIB{

	RBroydenFamily::RBroydenFamily(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		Initialization(prob, initialx, initialH, insoln);
	};

	void RBroydenFamily::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		SetProbX(prob, initialx, initialH, insoln);
		SetDefaultParams();
	};

	void RBroydenFamily::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		SolversLS::SetProbX(prob, initialx, insoln);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		bool initHisnull = (initialH == nullptr);
		if (initHisnull)
		{
			if (prob->GetDomain()->GetIsIntrinsic())
			{
				initialH = new LinearOPE(prob->GetDomain()->GetEMPTYINTR()->Getlength());
			}
			else
			{
				initialH = new LinearOPE(prob->GetDomain()->GetEMPTYEXTR()->Getlength());
			}
			initialH->ScaledIdOPE();
		}
		H = initialH->ConstructEmpty();
		tildeH = initialH->ConstructEmpty();
		initialH->CopyTo(H);
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();
		u = EMPTYETA->ConstructEmpty();
		if (initHisnull)
			delete initialH;
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RBroydenFamily::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		isconvex = false;
		nu = 1e-4;
		mu = 1;
		InitSteptype = QUADINTMOD;
		SolverName.assign("RBroydenFamily");
	};

	RBroydenFamily::~RBroydenFamily(void)
	{
		delete s;
		delete y;
		delete u;
		delete H;
		delete tildeH;
	};

	void RBroydenFamily::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RBROYDENFAMILY METHOD PARAMETERS:\n");
		status = (nu >= 0 && nu < 1) ? YES : NO;
		printf("nu            :%15g[%s],\t", nu, status);
		status = (mu >= 0) ? YES : NO;
		printf("mu            :%15g[%s]\n", mu, status);
		status = YES;
		printf("isconvex      :%15d[%s]\n", isconvex, status);
	};

	void RBroydenFamily::GetSearchDir(void)
	{
		HvRBroydenFamily(gf1, eta1);
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
	};

	void RBroydenFamily::UpdateData(void)
	{
		UpdateDataRBroydenFamily();
	};

	void RBroydenFamily::PrintInfo(void)
	{
		printf("\n\tbetay:%.3e,Phic:%.3e,inpss:%.3e,inpsy:%.3e,IsUpdateHessian:%d,", betay, phic, inpss, inpsy, isupdated);
		printf("\n");
	};
}; /*end of ROPTLIB namespace*/
