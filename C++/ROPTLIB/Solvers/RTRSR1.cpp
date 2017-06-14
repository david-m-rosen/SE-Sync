
#include "Solvers/RTRSR1.h"

/*Define the namespace*/
namespace ROPTLIB{

	RTRSR1::RTRSR1(const Problem *prob, const Variable *initialx, LinearOPE *initialB, const Variable *insoln)
	{
		Initialization(prob, initialx, initialB, insoln);
	};

	void RTRSR1::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialB, const Variable *insoln)
	{
		SetProbX(prob, initialx, initialB, insoln);
		SetDefaultParams();
	};

	void RTRSR1::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialB, const Variable *insoln)
	{
		SolversTR::SetProbX(prob, initialx, insoln);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		bool initBisnull = (initialB == nullptr);
		if (initBisnull)
		{
			if (prob->GetDomain()->GetIsIntrinsic())
			{
				initialB = new LinearOPE(prob->GetDomain()->GetEMPTYINTR()->Getlength());
			}
			else
			{
				initialB = new LinearOPE(prob->GetDomain()->GetEMPTYEXTR()->Getlength());
			}
			initialB->ScaledIdOPE();
		}
		B = initialB->ConstructEmpty();
		tildeB = initialB->ConstructEmpty();
		initialB->CopyTo(B);
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();
		prob->SetUseGrad(true);
		prob->SetUseHess(false);

		if (initBisnull)
			delete initialB;
	};

	void RTRSR1::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		theta = 0.1;
		kappa = 0.1;
		SolverName.assign("RTRSR1");
		isconvex = false;
	};

	RTRSR1::~RTRSR1(void)
	{
		delete s;
		delete y;
		delete B;
		delete tildeB;
	};

	void RTRSR1::HessianEta(Vector *Eta, Vector *result)
	{
		HvRTRSR1(Eta, result);
	};

	void RTRSR1::UpdateData(void)
	{
		UpdateDataRTRSR1();
	};

	void RTRSR1::Acceptence(void)
	{
		Mani->TranHInvTran(x1, eta2, x2, B, tildeB);
		tildeB->CopyTo(B);
	};

	void RTRSR1::CheckParams(void)
	{
		SolversTR::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RTRSR1 METHOD PARAMETERS:\n");
		status = YES;
		printf("isconvex      :%15d[%s]\n", isconvex, status);
	};

	void RTRSR1::PrintInfo(void)
	{
		printf("\n\tinpss:%.3e,IsUpdateHessian:%d,", inpss, isupdated);
		printf("\n");
	};
}; /*end of ROPTLIB namespace*/
