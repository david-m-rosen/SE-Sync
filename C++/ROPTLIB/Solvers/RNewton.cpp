#include "Solvers/RNewton.h"

/*Define the namespace*/
namespace ROPTLIB{

	RNewton::RNewton(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Initialization(prob, initialx, insoln);
	};

	void RNewton::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversLS::SetProbX(prob, initialx, insoln);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();

		r = EMPTYETA->ConstructEmpty();
		z = EMPTYETA->ConstructEmpty();
		delta = EMPTYETA->ConstructEmpty();
		Hd = EMPTYETA->ConstructEmpty();

		prob->SetUseGrad(true);
		prob->SetUseHess(true);
	};

	void RNewton::SetDefaultParams()
	{
		SolversLS::SetDefaultParams();
		theta = 1;
		kappa = 0.1;
		Min_Inner_Iter = 0;
		Max_Inner_Iter = 1000;
		useRand = false;
		InitSteptype = QUADINTMOD;

		SolverName.assign("RNewton");

		tCGLSstatusSetnames = new std::string[TCGLSSTATUSSETLENGTH];
		tCGLSstatusSetnames[LS_NEGCURVTURE].assign("NEGCURVTURE");
		tCGLSstatusSetnames[LS_LCON].assign("LCON");
		tCGLSstatusSetnames[LS_SCON].assign("SCON");
		tCGLSstatusSetnames[LS_MAXITER].assign("MAXITER");
	};

	void RNewton::CheckParams(void)
	{
		SolversLS::CheckParams();

		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("RNEWTON METHOD PARAMETERS:\n");
		status = (Min_Inner_Iter >= 0 && Min_Inner_Iter <= Max_Inner_Iter) ? YES : NO;
		printf("Min_Inner_Iter:%15d[%s],\t", Min_Inner_Iter, status);
		status = (Max_Inner_Iter >= 0 && Max_Inner_Iter >= Min_Inner_Iter) ? YES : NO;
		printf("Max_Inner_Iter:%15d[%s],\t", Max_Inner_Iter, status);
		status = (theta >= 1) ? YES : NO;
		printf("theta         :%15g[%s],\t", theta, status);
		status = (kappa > 0 && kappa < 1) ? YES : NO;
		printf("kappa         :%15g[%s]\n", kappa, status);
		printf("useRand       :%15d[%s]\n", useRand, status);
	};

	RNewton::~RNewton(void)
	{
		delete r;
		delete z;
		delete delta;
		delete Hd;
		delete[] tCGLSstatusSetnames;
	};

	void RNewton::GetSearchDir(void)
	{
		Mani->ScaleTimesVector(x1, 0, gf1, eta1);
		tCG_LS();
	};

	void RNewton::PreConditioner(Variable *x, Vector *eta, Vector *result)
	{
		// default one means no preconditioner.
		eta->CopyTo(result);
	};

	void RNewton::HessianEta(Vector *Eta, Vector *result) const
	{
		Prob->HessianEta(x1, Eta, result);
	};

	void RNewton::PrintInfo(void)
	{
		printf("\n\tnH:%d,tCGstatus:%s,innerIter:%d,", nH, tCGLSstatusSetnames[tCGLSstatus].c_str(), innerIter);
		printf("\n");
	};

	void RNewton::tCG_LS(void)
	{
		double e_Pe, r_r, norm_r, norm_r0, d_Pd, z_r, e_Pd, d_Hd, alphatemp, e_Pe_new, tempnum, zold_rold, betatemp;
		integer j;

		if (useRand)
		{
			HessianEta(eta1, r);
			Mani->VectorAddVector(x1, gf1, r, r);
			e_Pe = Mani->Metric(x1, eta1, eta1);
		}
		else
		{
			gf1->CopyTo(r);
			e_Pe = 0;
		}

		r_r = Mani->Metric(x1, r, r);
		norm_r = sqrt(r_r);
		norm_r0 = norm_r;

		PreConditioner(x1, r, z);

		z_r = Mani->Metric(x1, z, r);
		d_Pd = z_r;

		Mani->ScaleTimesVector(x1, -1.0, z, delta);

		if (useRand)
			e_Pd = Mani->Metric(x1, eta1, delta);
		else
			e_Pd = 0;

		tCGLSstatus = LS_MAXITER; /* pre-assume termination j == max_inner*/
		for (j = 0; j < Max_Inner_Iter; j++)
		{
			HessianEta(delta, Hd); nH++;
			d_Hd = Mani->Metric(x1, delta, Hd);
			alphatemp = z_r / d_Hd;
			e_Pe_new = e_Pe + 2.0 * alphatemp * e_Pd + alphatemp * alphatemp * d_Pd;

			if (d_Hd <= 0)
			{
				Mani->VectorAddVector(x1, delta, eta1, eta1);
				tCGLSstatus = LS_NEGCURVTURE; /* negative curvature*/
				break;
			}
			e_Pe = e_Pe_new;
			Mani->scalarVectorAddVector(x1, alphatemp, delta, eta1, eta1);

			Mani->scalarVectorAddVector(x1, alphatemp, Hd, r, r);

			Mani->Projection(x1, r, zeta);
			zeta->CopyTo(r);

			r_r = Mani->Metric(x1, r, r);
			norm_r = sqrt(r_r);

			tempnum = pow(norm_r0, theta);
			if (j >= Min_Inner_Iter && norm_r <= norm_r0 * ((tempnum < kappa) ? tempnum : kappa))
			{
				if (kappa < tempnum)
					tCGLSstatus = LS_LCON; /* linear convergence*/
				else
					tCGLSstatus = LS_SCON; /* superlinear convergence*/
				break;
			}

			PreConditioner(x1, r, z);

			zold_rold = z_r;
			z_r = Mani->Metric(x1, z, r);

			betatemp = z_r / zold_rold;
			Mani->scalarVectorMinusVector(x1, betatemp, delta, z, delta);
			e_Pd = betatemp * (e_Pd + alphatemp * d_Pd);
			d_Pd = z_r + betatemp * betatemp * d_Pd;
		}
		innerIter = j;
	};

	void RNewton::SetParams(PARAMSMAP params)
	{
		SolversLS::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("useRand"))
			{
				useRand = ((static_cast<integer> (iter->second)) != 0);
			}
			else
				if (iter->first == static_cast<std::string> ("Max_Inner_Iter"))
				{
					Max_Inner_Iter = static_cast<integer> (iter->second);
				}
				else
					if (iter->first == static_cast<std::string> ("Min_Inner_Iter"))
					{
						Min_Inner_Iter = static_cast<integer> (iter->second);
					}
					else
						if (iter->first == static_cast<std::string> ("theta"))
						{
							theta = iter->second;
						}
						else
							if (iter->first == static_cast<std::string> ("kappa"))
							{
								kappa = static_cast<integer> (iter->second);
							}
		}
	};
}; /*end of ROPTLIB namespace*/
