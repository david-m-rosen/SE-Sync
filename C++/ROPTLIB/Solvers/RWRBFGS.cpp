
#include "Solvers/RWRBFGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	RWRBFGS::RWRBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		Initialization(prob, initialx, initialH, insoln);
	};

	void RWRBFGS::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		SetProbX(prob, initialx, initialH, insoln);
		SetDefaultParams();
	};

	void RWRBFGS::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
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

		prob->SetUseGrad(true);
		prob->SetUseHess(false);
		if (initHisnull)
			delete initialH;
	};

	void RWRBFGS::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		isconvex = false;
		nu = 1e-4;
		mu = 1;
		InitSteptype = QUADINTMOD;
		SolverName.assign("RWRBFGS");
	};

	RWRBFGS::~RWRBFGS(void)
	{
		delete s;
		delete y;
		delete H;
		delete tildeH;
	};

	void RWRBFGS::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RWRBFGS METHOD PARAMETERS:\n");
		status = (nu >= 0 && nu < 1) ? YES : NO;
		printf("nu            :%15g[%s],\t", nu, status);
		status = (mu >= 0) ? YES : NO;
		printf("mu            :%15g[%s]\n", mu, status);
		status = YES;
		printf("isconvex      :%15d[%s]\n", isconvex, status);
	};

	void RWRBFGS::GetSearchDir(void)
	{
		HvRWRBFGS(gf1, eta1);
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
	};

	void RWRBFGS::UpdateData(void)
	{
		UpdateDataRWRBFGS();
	};

	void RWRBFGS::PrintInfo(void)
	{
		printf("\n\tinpss:%.3e,inpsy:%.3e,IsUpdateHessian:%d,", inpss, inpsy, isupdated);
		printf("\n");
	};
}; /*end of ROPTLIB namespace*/
