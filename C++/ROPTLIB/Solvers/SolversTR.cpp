
#include "Solvers/SolversTR.h"

/*Define the namespace*/
namespace ROPTLIB{

	void SolversTR::Run(void)
	{
		Variable *xTemp;
		Vector *gfTemp;
		Solvers::Run();
		starttime = getTickCount();
		double sqeps = sqrt(std::numeric_limits<double>::epsilon());

		f1 = Prob->f(x1); nf++;
		f2 = f1;
		Prob->Grad(x1, gf1); ng++;

		ngf0 = sqrt(Mani->Metric(x1, gf1, gf1));
		ngf = ngf0;

		iter = 0;
		if (Debug >= ITERRESULT)
		{
			printf("i:%d,f:%.3e,|gf|:%.3e,\n", iter, f1, ngf);
			timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS;
			funSeries[iter] = f1; gradSeries[iter] = ngf;
			distSeries[iter] = ((soln == nullptr) ? 0 : Mani->Dist(x1, soln));
		}
		Delta = initial_Delta;
		bool isstop = IsStopped();
		while (((! isstop) && iter < Max_Iteration) || iter < Min_Iteration)
		{
			InitialVector(); // Obtain initial guess, eta1, for local model
			tCG_TR(); // obtain eta2
			Mani->Retraction(x1, eta2, x2);	nR++;
			f2 = Prob->f(x2); nf++;
			HessianEta(eta2, zeta); nH++; // Hessian * eta2

			Mani->scalarVectorAddVector(x1, 0.5, zeta, gf1, eta1);
			rho = (f1 - f2) / (-Mani->Metric(x1, eta2, eta1));
			UpdateData(); // Update S&Y or H or B

			if (rho > 0.75)
			{
				if (tCGstatus == TR_EXCREGION || tCGstatus == TR_NEGCURVTURE)
					Delta *= Magnified_tau;
				if (Delta > maximum_Delta)
				{
					if (Debug >= FINALRESULT)
					{
						printf("reach the maximum of radius\n");
					}
					Delta = maximum_Delta;
				}
			}
			else
			if (rho < 0.25)
			{
				Delta *= Shrinked_tau;
				if (Delta < minimum_Delta)
				{
					if (Debug >= FINALRESULT)
					{
						printf("reach the minimum of radius\n");
					}
					break;
				}
			}
			if (rho > Acceptence_Rho || (fabs(f1 - f2) / (fabs(f1) + 1) < sqeps && f2 < f1))
			{
				Acceptence(); // Algorithm specific operations
				ngf = sqrt(Mani->Metric(x2, gf2, gf2));
				isstop = IsStopped(); /*This is done when the candidate is accepted. This is necessary for partly smooth stopping criterion*/
				xTemp = x1; x1 = x2; x2 = xTemp;
				gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;
				iter++;
				if (Debug >= ITERRESULT && iter % OutputGap == 0)
				{
					PrintGenInfo();
					PrintInfo(); // Output information specific to Algorithms
				}
				f1 = f2;
			}
			else
			{
				iter++;
				if (Debug >= ITERRESULT && iter % OutputGap == 0)
				{
					printf("X_{%d} WAS REJECTED.\n", iter);
					PrintGenInfo();
					PrintInfo(); // Output information specific to Algorithms
				}
			}

			if (Debug >= ITERRESULT)
			{
				timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS;
				funSeries[iter] = f2; gradSeries[iter] = ngf;
				distSeries[iter] = ((soln == nullptr) ? 0 : Mani->Dist(x1, soln));
			}
		}
		ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;
		if (Debug >= ITERRESULT)
			lengthSeries = iter + 1;

		if (Debug >= FINALRESULT)
		{
			printf("Iter:%d,f:%.3e,", iter, f2);
			if (nsubgf != -1)
			{
				printf("nsubgf:%.3e,", nsubgf);
			}
			printf("|gf|:%.3e,|gf|/|gf0|:%.3e,time:%.2e,nf:%d,ng:%d,nR:%d,", ngf, ngf / ngf0, ComTime, nf, ng, nR);
			if (nH != 0)
			{
				printf("nH:%d,", nH);
			}
			if (nV != 0)
			{
				printf("nV(nVp):%d(%d),", nV, nVp);
			}
			printf("\n");
		}
	};

	void SolversTR::InitialVector(void)
	{
		Mani->ScaleTimesVector(x1, 0, gf1, eta1);
	};

	void SolversTR::tCG_TR(void)
	{
		double e_Pe, r_r, norm_r, norm_r0, d_Pd, z_r, e_Pd, d_Hd, alphatemp, e_Pe_new, tempnum, zold_rold, betatemp, tautemp;
		integer j;

		if (useRand)
		{
			HessianEta(eta1, r); nH++;
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

		tCGstatus = TR_MAXITER; /* pre-assume termination j == max_inner*/

		eta1->CopyTo(eta2);

		for (j = 0; j < Max_Inner_Iter; j++)
		{
			HessianEta(delta, Hd); nH++;
			d_Hd = Mani->Metric(x1, delta, Hd);
			alphatemp = z_r / d_Hd;
			e_Pe_new = e_Pe + 2.0 * alphatemp * e_Pd + alphatemp * alphatemp * d_Pd;

			if (d_Hd <= 0 || e_Pe_new >= (Delta * Delta))
			{
				tautemp = (-e_Pd + sqrt(e_Pd * e_Pd + d_Pd * (Delta * Delta - e_Pe))) / d_Pd;
				Mani->scalarVectorAddVector(x1, tautemp, delta, eta2, eta2);

				if (d_Hd < 0)
					tCGstatus = TR_NEGCURVTURE; /* negative curvature*/
				else
					tCGstatus = TR_EXCREGION; /* exceeded trust region*/
				break;
			}
			e_Pe = e_Pe_new;
			Mani->scalarVectorAddVector(x1, alphatemp, delta, eta2, eta2);

			Mani->scalarVectorAddVector(x1, alphatemp, Hd, r, r);

			Mani->Projection(x1, r, zeta);
			zeta->CopyTo(r);

			r_r = Mani->Metric(x1, r, r);
			norm_r = sqrt(r_r);

			tempnum = pow(norm_r0, theta);
			if (j >= Min_Inner_Iter && norm_r <= norm_r0 * ((tempnum < kappa) ? tempnum : kappa))
			{
				if (kappa < tempnum)
					tCGstatus = TR_LCON; /* linear convergence*/
				else
					tCGstatus = TR_SCON; /* superlinear convergence*/
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

	void SolversTR::PreConditioner(Variable *x, Vector *eta, Vector *result)
	{
		// default one means no preconditioner.
		eta->CopyTo(result);
	};

	void SolversTR::PrintGenInfo(void)
	{
		Solvers::PrintGenInfo();
		printf("nH:%d,rho:%.2e,radius:%.3e,tCGstatus:%s,innerIter:%d,", nH, rho, Delta, tCGstatusSetnames[tCGstatus].c_str(), innerIter);
	};

	void SolversTR::CheckParams(void)
	{
		Solvers::CheckParams();

		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("TRUST REGION TYPE METHODS PARAMETERS:\n");
		status = (initial_Delta > 0) ? YES : NO;
		printf("initial_Delta :%15g[%s],\t", initial_Delta, status);
		status = (Acceptence_Rho > 0 && Acceptence_Rho < 0.25) ? YES : NO;
		printf("Acceptence_Rho:%15g[%s]\n", Acceptence_Rho, status);
		status = (Shrinked_tau > 0 && Shrinked_tau < 1) ? YES : NO;
		printf("Shrinked_tau  :%15g[%s],\t", Shrinked_tau, status);
		status = (Magnified_tau > 1) ? YES : NO;
		printf("Magnified tau :%15g[%s]\n", Magnified_tau, status);
		status = (minimum_Delta > 0 && minimum_Delta <= maximum_Delta) ? YES : NO;
		printf("minimum_Delta :%15g[%s],\t", minimum_Delta, status);
		status = (maximum_Delta > 0 && maximum_Delta >= minimum_Delta) ? YES : NO;
		printf("maximum_Delta :%15g[%s]\n", maximum_Delta, status);
		status = (Min_Inner_Iter >= 0 && Min_Inner_Iter <= Max_Inner_Iter) ? YES : NO;
		printf("Min_Inner_Iter:%15d[%s],\t", Min_Inner_Iter, status);
		status = (Max_Inner_Iter >= 0 && Max_Inner_Iter >= Min_Inner_Iter) ? YES : NO;
		printf("Max_Inner_Iter:%15d[%s],\t", Max_Inner_Iter, status);
		status = (theta >= 0) ? YES : NO;
		printf("theta         :%15g[%s],\t", theta, status);
		status = (kappa > 0 && kappa < 1) ? YES : NO;
		printf("kappa         :%15g[%s]\n", kappa, status);
		status = YES;
		printf("useRand       :%15d[%s]\n", useRand, status);
	};

	void SolversTR::UpdateData(void)
	{
	};

	void SolversTR::Acceptence(void)
	{
		Prob->Grad(x2, gf2); ng++;
	};

	void SolversTR::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		Solvers::SetProbX(prob, initialx, insoln);
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();

		r = EMPTYETA->ConstructEmpty();
		z = EMPTYETA->ConstructEmpty();
		delta = EMPTYETA->ConstructEmpty();
		Hd = EMPTYETA->ConstructEmpty();
	};

	void SolversTR::SetDefaultParams()
	{
		Solvers::SetDefaultParams();
		nH = 0;
		Acceptence_Rho = 0.1;
		Shrinked_tau = 0.25;
		Magnified_tau = 2;
		minimum_Delta = std::numeric_limits<double>::epsilon();
		maximum_Delta = 10000;
		useRand = false;
		Max_Inner_Iter = 1000;
		Min_Inner_Iter = 0;
		theta = 1;
		kappa = 0.1;
		initial_Delta = 1;
		tCGstatusSetnames = new std::string[TCGSTATUSSETLENGTH];
		tCGstatusSetnames[TR_NEGCURVTURE].assign("NEGCURVTURE");
		tCGstatusSetnames[TR_EXCREGION].assign("EXCREGION");
		tCGstatusSetnames[TR_LCON].assign("LCON");
		tCGstatusSetnames[TR_SCON].assign("SCON");
		tCGstatusSetnames[TR_MAXITER].assign("MAXITER");
	};

	SolversTR::~SolversTR(void)
	{
		delete r;
		delete z;
		delete delta;
		delete Hd;

		delete[] tCGstatusSetnames;
	};

	void SolversTR::SetParams(PARAMSMAP params)
	{
		QuasiNewton::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Acceptence_Rho"))
			{
				Acceptence_Rho = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Shrinked_tau"))
			{
				Shrinked_tau = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Magnified_tau"))
			{
				Magnified_tau = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("minimum_Delta"))
			{
				minimum_Delta = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("maximum_Delta"))
			{
				maximum_Delta = iter->second;
			}
			else
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
				kappa = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("initial_Delta"))
			{
				initial_Delta = iter->second;
			}
		}
	};
}; /*end of ROPTLIB namespace*/
