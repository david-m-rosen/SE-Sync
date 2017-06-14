
#include "Solvers/RBFGS.h"

/*Define the namespace*/
namespace ROPTLIB{

	RBFGS::RBFGS(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		Initialization(prob, initialx, initialH, insoln);
	};

	void RBFGS::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		SetProbX(prob, initialx, initialH, insoln);
		SetDefaultParams();
	};

	void RBFGS::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
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

		if (initHisnull)
			delete initialH;
		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void RBFGS::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		isconvex = false;
		nu = 1e-4;
		mu = 1;
		InitSteptype = QUADINTMOD;
		SolversLS::SolverName.assign("RBFGS");
	};

	RBFGS::~RBFGS(void)
	{
		delete s;
		delete y;
		delete H;
		delete tildeH;
	};

	void RBFGS::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RBFGS METHOD PARAMETERS:\n");
		status = (nu >= 0 && nu < 1) ? YES : NO;
		printf("nu            :%15g[%s],\t", nu, status);
		status = (mu >= 0) ? YES : NO;
		printf("mu            :%15g[%s]\n", mu, status);
		status = YES;
		printf("isconvex      :%15d[%s]\n", isconvex, status);
	};

	void RBFGS::GetSearchDir(void)
	{
		HvRBFGS(gf1, eta1);
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
	};

	void RBFGS::UpdateData(void)
	{
		UpdateDataRBFGS();
	};

	void RBFGS::PrintInfo(void)
	{
		printf("\n\tbetay:%.3e,inpss:%.3e,inpsy:%.3e,IsUpdateHessian:%d,", betay, inpss, inpsy, isupdated);
		printf("\n");
	};

}; /*end of ROPTLIB namespace*/
