
#include "Solvers/RBFGSLPSub.h"

/*Define the namespace*/
namespace ROPTLIB{

	RBFGSLPSub::RBFGSLPSub(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		Initialization(prob, initialx, initialH, insoln);
	};

	void RBFGSLPSub::Initialization(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		SetProbX(prob, initialx, initialH, insoln);
		SetDefaultParams();
	};

	void RBFGSLPSub::SetProbX(const Problem *prob, const Variable *initialx, LinearOPE *initialH, const Variable *insoln)
	{
		SolversLSLPSub::SetProbX(prob, initialx, insoln);

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

	void RBFGSLPSub::SetDefaultParams()
	{
		SolversLSLPSub::SetDefaultParams();
		isconvex = false;
		LineSearch_LS = WOLFELP;
		InitSteptype = ONESTEP;
		lambdaLower = 1e-2;//-- 7;
		lambdaUpper = 1e2;//-- 7;
		Hv = &QuasiNewton::HvRBFGSSub;
		SolversLSLPSub::SolverName.assign("RBFGSLPSub");
	};

	RBFGSLPSub::~RBFGSLPSub(void)
	{
		delete s;
		delete y;
		delete H;
		delete tildeH;
	};

	void RBFGSLPSub::PrintInfo(void)
	{
		printf("\n\tbetay:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", betay, inpss, inpsy, inpyy, isupdated);
		printf("\n");
	};

	void RBFGSLPSub::CheckParams(void)
	{
		SolversLSLPSub::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RBFGSLPSub METHOD PARAMETERS:\n");
		status = (lambdaLower > 0 && lambdaLower < lambdaUpper) ? YES : NO;
		printf("lambdaLower   :%15g[%s],\t", lambdaLower, status);
		status = (lambdaUpper >= lambdaLower) ? YES : NO;
		printf("lambdaUpper   :%15g[%s]\n", lambdaUpper, status);
		status = YES;
		printf("isconvex      :%15d[%s]\n", isconvex, status);
	};

	void RBFGSLPSub::UpdateData(void)
	{
		gf->CopyTo(gf1);
		UpdateDataRBFGSSub();
	};
}; /*end of ROPTLIB namespace*/
