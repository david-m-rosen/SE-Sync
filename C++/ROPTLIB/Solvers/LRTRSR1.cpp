
#include "Solvers/LRTRSR1.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRTRSR1::LRTRSR1(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void LRTRSR1::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversTR::SetProbX(prob, initialx, insoln);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();
		s = EMPTYETA->ConstructEmpty();
		y = EMPTYETA->ConstructEmpty();

		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void LRTRSR1::SetDefaultParams()
	{
		SolversTR::SetDefaultParams();
		theta = 0.1;
		kappa = 0.1;
		isconvex = false;
		LengthSY = 4;
		S = nullptr;
		Y = nullptr;
		YMGS = nullptr;
		inpss = 0;
		inpsy = 0;
		inpyy = 0;
		Currentlength = 0;
		beginidx = 0;
		SS = nullptr;
		SY = nullptr;
		PMGQ = nullptr;
		P = nullptr;
		gamma = 1;
		ischangedSandY = true;
		SolverName.assign("LRTRSR1");
	};

	LRTRSR1::~LRTRSR1(void)
	{
		delete s;
		delete y;
		DeleteVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		if (SS != nullptr)
			delete[] SS;
		if (SY != nullptr)
			delete[] SY;
		if (PMGQ != nullptr)
			delete[] PMGQ;
		if (P != nullptr)
			delete[] P;
	};

	void LRTRSR1::Run(void)
	{
		DeleteVectors(S, LengthSY);
		NewVectors(S, LengthSY);
		DeleteVectors(Y, LengthSY);
		NewVectors(Y, LengthSY);
		DeleteVectors(YMGS, LengthSY);
		NewVectors(YMGS, LengthSY);
		if (SS != nullptr)
			delete[] SS;
		SS = new double[LengthSY * LengthSY];
		if (SY != nullptr)
			delete[] SY;
		SY = new double[LengthSY * LengthSY];
		if (PMGQ != nullptr)
			delete[] PMGQ;
		PMGQ = new double[LengthSY * LengthSY];
		if (P != nullptr)
			delete[] P;
		P = new integer[LengthSY];
		SolversTR::Run();
	};

	void LRTRSR1::CheckParams(void)
	{
		SolversTR::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("LRTRSR1 METHOD PARAMETERS:\n");
		status = YES;
		printf("isconvex      :%15d[%s],\t", isconvex, status);
		status = (LengthSY >= 0) ? YES : NO;
		printf("LengthSY      :%15d[%s]\n", LengthSY, status);
	};

	void LRTRSR1::HessianEta(Vector *Eta, Vector *result)
	{
		HvLRTRSR1(Eta, result);
	};

	void LRTRSR1::UpdateData(void)
	{
		UpdateDataLRTRSR1();
	};

	void LRTRSR1::Acceptence(void)
	{
		for (integer i = 0; i < Currentlength; i++)
		{
			Mani->VectorTransport(x1, eta2, x2, S[i], S[i]);
			Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]);
		}
		ischangedSandY = true;
	};

	void LRTRSR1::PrintInfo(void)
	{
		printf("\n\tgamma:%.3e,inpss:%.3e,inpsy:%.3e,inpyy:%.3e,IsUpdateHessian:%d,", gamma, inpss, inpsy, inpyy, isupdated);
		printf("\n");
	};
}; /*end of ROPTLIB namespace*/
