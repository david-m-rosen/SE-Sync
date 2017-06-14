
#include "Solvers/SolversLSLPSub.h"

/*Define the namespace*/
namespace ROPTLIB{
	void SolversLSLPSub::Run(void)
	{
		Variable *xTemp;
		Vector *gfTemp;
		pre_funs.clear();
		Solvers::Run();
		ChooseLinesearch();

		LSstatus = SUCCESS;
		f1 = Prob->f(x1); nf++;
		f2 = f1;
		Prob->Grad(x1, gf1); ng++;
		ngf0 = sqrt(Mani->Metric(x1, gf1, gf1));
		ngf = ngf0;
		newslope = 0;
		iter = 0;
		if (Debug >= ITERRESULT)
		{
			printf("i:%d,f:%.3e,\n", iter, f1);
			timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS;
			funSeries[iter] = f1;
			gradSeries[iter] = ngf;
			distSeries[iter] = ((soln == nullptr) ? 0 : Mani->Dist(x1, soln));
		}
		bool isstop = false;

		/*If the intrinsic representation is used, then exeta1 and exeta2 are used to store the extrinsic representations of eta1 and eta2 respectively*/
		if (Prob->GetDomain()->GetIsIntrinsic())
		{
			exeta1 = Prob->GetDomain()->GetEMPTYEXTR()->ConstructEmpty();
			exeta2 = Prob->GetDomain()->GetEMPTYEXTR()->ConstructEmpty();
		}

		/*Start the loop*/
		while ((((!isstop) && iter < Max_Iteration) || iter < Min_Iteration) && LSstatus == SUCCESS)
		{
			GetSearchDir(); // Obtain search direction eta1 and gf

			/*Call the function to check whether the stopping criterion is satisfied or not.
			The default function is written in Solvers.h and Solvers.cpp*/
			isstop = IsStopped();

			if (fabs(ngf) <= Del || fabs(ngf) <= Tolerance)
			{
				Eps *= Theta_eps;
				Eps = (Eps > Min_Eps) ? Eps : Min_Eps;
				if (Eps > Min_Eps)
				{
					Del *= Theta_del;
				}
				else
				{
					while (fabs(ngf) < Del)
					{
						Del *= Theta_del;
					}
				}
				if (Debug >= ITERRESULT)
					printf("Shinking Epsilon and Delta to %g and %g respectively, |gf|:%g.\n", Eps, Del, ngf);
				continue;
			}

			initialslope = Mani->Metric(x1, gf, eta1);
			/*Compute initial step size for the next iteration*/
			InitialStepSize();

			initiallength = stepsize;

			/*Start a line search algorithm. If the intrinsic representation is used for the search direction, then converting it into the extrinsic representation.*/
			if (Prob->GetDomain()->GetIsIntrinsic())
			{
				Prob->GetDomain()->ObtainExtr(x1, eta1, exeta1);
			}

			/* Call the specified linesearch algorithm.
			Note that in the linesearch algorithm, we need to obtain
			accepted stepsize, eta2=stepsize*eta1, x2 = R_{x_1}(eta_2), f2 = f(x2), and gf2 = grad f(x_2) */
			if (LineSearch_LS == INPUTFUN)
			{
				if (Prob->GetDomain()->GetIsIntrinsic())
				{
					//Vector *exeta1 = Prob->GetDomain()->GetEMPTYEXTR()->ConstructEmpty();
					//Prob->GetDomain()->ObtainExtr(x1, eta1, exeta1);
					stepsize = LinesearchInput(iter, x1, exeta1, stepsize, initialslope, Prob, this);
					//delete exeta1;
				}
				else
				{
					stepsize = LinesearchInput(iter, x1, eta1, stepsize, initialslope, Prob, this);
				}
				stepsize = (stepsize < std::numeric_limits<double>::epsilon()) ? initiallength : stepsize;
				initiallength = stepsize;
				if (!IsPureLSInput)
				{
					SolversLS::LinesearchArmijo();
				}
				else
				{
					f2 = h(); nf++;
					Prob->Grad(x2, gf2); ng++;
				}
			}
			else
			{
				(this->*Linesearch)();
			}

			/*Output debug information if necessary.*/
			if (LSstatus < SUCCESS && Debug >= FINALRESULT)
			{
				printf("Linesearch fails! LSstatus:%s\n", LSstatusSetnames[LSstatus].c_str());
			}

			iter++;

			/*Update the Hessian approximation for quasi-Newton methods
			or obtain search direction candadite for Riemannian nonlinear conjugate gradient*/
			UpdateData();

			if (Debug >= ITERRESULT)
			{
				/*Output information*/
				if (iter % OutputGap == 0)
				{
					PrintGenInfo();
					PrintInfo(); // Output information specific to Algorithms
				}
				/*Store debug information in the arrays*/
				timeSeries[iter] = static_cast<double>(getTickCount() - starttime) / CLK_PS;
				funSeries[iter] = f2; gradSeries[iter] = ngf;
				distSeries[iter] = ((soln == nullptr) ? 0 : Mani->Dist(x2, soln));
			}

			/*Switch information at x1 and x2*/
			xTemp = x1; x1 = x2; x2 = xTemp;
			gfTemp = gf1; gf1 = gf2; gf2 = gfTemp;
			pre_funs.push_front(f1);
			if (pre_funs.size() > Num_pre_funs && pre_funs.size() > 1)
				pre_funs.pop_back();
			f1 = f2;
		}
		if (Prob->GetDomain()->GetIsIntrinsic())
		{
			delete exeta1;
			delete exeta2;
		}
		ComTime = static_cast<double>(getTickCount() - starttime) / CLK_PS;
		if (Debug >= ITERRESULT)
			lengthSeries = iter + 1;
		if (Debug >= FINALRESULT)
		{
			printf("Iter:%d,f:%.3e,|gf|:%.3e,time:%.2e,nsubprob:%d,nf:%d,ng:%d,nR:%d,", iter, f2, ngf, ComTime, subprobtimes, nf, ng, nR);
			if (nH != 0)
			{
				printf("nH:%d,", nH);
			}
			if (nV != 0)
			{
				printf("nV(nVp):%d(%d),", nV, nVp);
			}
			printf("Eps:%.3e,", Eps);
			printf("\n");
		}
	};

	bool SolversLSLPSub::IsStopped(void)
	{
		if (ngf < Tolerance && Eps <= Min_Eps + std::numeric_limits<double>::epsilon() * 10)// (Eps <= ((Tolerance < Min_Eps) ? Min_Eps : Tolerance)))
			return true;
		return false;
	};

	void SolversLSLPSub::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SolversLS::SetProbX(prob, initialx, insoln);

		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();

		gf = EMPTYETA->ConstructEmpty();

		prob->SetUseGrad(true);
		prob->SetUseHess(false);
	};

	void SolversLSLPSub::SetDefaultParams(void)
	{
		SolversLS::SetDefaultParams();
		Theta_eps = 0.01;
		Theta_del = 0.01;
		Eps = 1;
		Min_Eps = 1e-6;
		Del = 1;

		Hv = nullptr;
	};

	SolversLSLPSub::~SolversLSLPSub(void)
	{
		delete gf;
	};

	void SolversLSLPSub::CheckParams(void)
	{
		SolversLS::CheckParams();
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;

		printf("SolversLSLPSub METHOD PARAMETERS:\n");
		status = (Eps > 0) ? YES : NO;
		printf("Eps           :%15g[%s],\t", Eps, status);
		status = (Theta_eps > 0 && Theta_eps < 1) ? YES : NO;
		printf("Theta_eps     :%15g[%s]\n", Theta_eps, status);
		status = (Min_Eps > 0 && Min_Eps < 1) ? YES : NO;
		printf("Min_Eps       :%15g[%s],\t", Min_Eps, status);
		status = (Del > 0) ? YES : NO;
		printf("Del           :%15g[%s]\n", Del, status);
		status = (Theta_del > 0 && Theta_del < 1) ? YES : NO;
		printf("Theta_del     :%15g[%s]\n", Theta_del, status);
	};

	void SolversLSLPSub::GetSearchDir(void)
	{
		gf1->CopyTo(gfs[0]);
		Currentlengthgfs = 1;
		Sub_soln = new double[Lengthgfs];
		for (integer i = 0; i < Lengthgfs; i++)
			Sub_soln[i] = 1;
		double Pngfsq = 0, neta1 = 0, hb;
		double *YtPY = new double[Lengthgfs * Lengthgfs];
		for (integer i = 0; i < Lengthgfs - 1; i++)
		{
			/*gf = argmin \| v \|_P^2, v in convex hall of W, let eta1 = - P gf.*/
			Pngfsq = MinPNormConHull(Mani, x1, gfs, Currentlengthgfs, this, Hv, gf, YtPY, Lengthgfs);
			(this->*Hv)(gf, eta1);
			Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
			if (Currentlengthgfs > 1)
				subprobtimes++;

			//Pngfsq = MinPnormWv();

			if (sqrt(Pngfsq) > Del * lambdaLower)
				ngf = sqrt(Mani->Metric(x1, gf, gf));
			else
			{
				ngf = sqrt(MinPNormConHull(Mani, x1, gfs, Currentlengthgfs, nullptr, nullptr, gf, nullptr, 0));
				subprobtimes++;
			}
			if (ngf <= Del)
			{
				delete[] Sub_soln;
				delete[] YtPY;
				return;
			}
			neta1 = sqrt(Mani->Metric(x1, eta1, eta1));
			stepsize = Eps / neta1;
			if (Prob->GetDomain()->GetIsIntrinsic())
			{
				Prob->GetDomain()->ObtainExtr(x1, eta1, exeta1);
			}
			f2 = h();
			hb = f2 - f1 + LS_alpha * Eps / neta1 * Pngfsq;
			if (hb <= 0)
			{
				delete[] Sub_soln;
				delete[] YtPY;
				return;
			}
			/*get gf2 and stepsize*/
			Increasing(neta1, Pngfsq, hb);
			Mani->InverseVectorTransport(x1, eta2, x2, gf2, gf2);
			betay = Mani->Beta(x1, eta2);
			Mani->ScaleTimesVector(x1, 1.0 / betay, gf2, gf2);
			gf2->CopyTo(gfs[Currentlengthgfs]);
			Currentlengthgfs++;
			if (Currentlengthgfs == Lengthgfs && Debug >= ITERRESULT)
			{
				printf("Warning: the number of W reaches its upper-bound: %d!\n", Currentlengthgfs);
			}
		}
		//if (Mani->Metric(x1, gf1, eta1) > -std::numeric_limits<double>::epsilon())
		//{
		//	Mani->ScaleTimesVector(x1, -1.0, gf1, eta1);
		//}
		delete[] Sub_soln;
		delete[] YtPY;
	};

	void SolversLSLPSub::Increasing(double neta1, double Pngfsq, double hb)
	{
		double a = 0, b = stepsize, ht;
		integer times = 0;
		while (times < 10)
		{
			if (dh() + LS_alpha * Pngfsq < 0)
			{
				stepsize = (a + b) / 2;
				f2 = h();
				ht = f2 - f1 + LS_alpha * Eps / neta1 * Pngfsq;
				if (hb > ht)
				{
					a = stepsize;
				}
				else
				{
					b = stepsize;
				}
			}
			else
			{
				break;
			}
			times++;
		}
		if (times == 10 && Debug >= ITERRESULT)
		{
			printf("warning: the loop in SolversLSLPSub::Increasing reaches the upperbound!\n");
		}
	};

	double SolversLSLPSub::MinPnormWv(void)
	{
		if (Currentlengthgfs == 1)
		{
			gfs[0]->CopyTo(gf);
			(this->*Hv)(gfs[0], eta1);
			const double *PWptr = eta1->ObtainReadData();
			const double *Wptr = gfs[0]->ObtainReadData();
			integer Length = eta1->Getlength();
			double result = ddot_(&Length, const_cast<double *> (Wptr), &GLOBAL::IONE, const_cast<double *> (PWptr), &GLOBAL::IONE);
			Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
			return result;
		}

		SphereConvexHull subprob(Mani, x1, gfs, Currentlengthgfs, this, Hv);
		Sphere Domain(Currentlengthgfs);
		subprob.SetDomain(&Domain);
		SphereVariable InitX(Currentlengthgfs);

		double *InitXptr = InitX.ObtainWriteEntireData();
		double deno = sqrt((static_cast<double>(Currentlengthgfs - 1)) / Currentlengthgfs);
		for (integer i = 0; i < InitX.Getlength() - 1; i++)
			InitXptr[i] = sqrt(Sub_soln[i]) * deno;
		InitXptr[InitX.Getlength() - 1] = 1.0 / sqrt(static_cast<double>(Currentlengthgfs));

		Solvers *RTRNewtonsolver = new RTRNewton(&subprob, &InitX);
		RTRNewtonsolver->Stop_Criterion = GRAD_F;
		RTRNewtonsolver->Debug = NOOUTPUT;//--- FINALRESULT;
		//RTRNewtonsolver->Max_Iteration = 20;
		RTRNewtonsolver->Tolerance = 1e-12;
		RTRNewtonsolver->Run(); 
		if (Currentlengthgfs > 1)
			subprobtimes++;
		double result = RTRNewtonsolver->Getfinalfun();
		const Variable *soln = RTRNewtonsolver->GetXopt();
		const double *solnptr = soln->ObtainReadData();
		for (integer i = 0; i < Currentlengthgfs; i++)
			Sub_soln[i] = solnptr[i] * solnptr[i];
		const SharedSpace *SharedWxsq = soln->ObtainReadTempData("Wxsq");
		SharedWxsq->GetSharedElement()->CopyTo(gf);
		const SharedSpace *SharedPWxsq = soln->ObtainReadTempData("PWxsq");
		SharedPWxsq->GetSharedElement()->CopyTo(eta1);
		Mani->ScaleTimesVector(x1, -1.0, eta1, eta1);
		delete RTRNewtonsolver;
		return result;
	};

	void SolversLSLPSub::SetParams(PARAMSMAP params)
	{
		SolversLS::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Eps"))
			{
				Eps = static_cast<double> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Theta_eps"))
			{
				Theta_eps = static_cast<double> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Min_Eps"))
			{
				Min_Eps = static_cast<double> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Del"))
			{
				Del = static_cast<double> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Theta_del"))
			{
				Theta_del = static_cast<double> (iter->second);
			}
		}
	};
}; /*end of ROPTLIB namespace*/
