
#include "Solvers/Solvers.h"

/*Define the namespace*/
namespace ROPTLIB{

	void Solvers::OutPutResults(Variable *inx1, double &inf1, double &inngf0, double &inngf, integer &initer,
		integer &innf, integer &inng, integer &innR, integer &innV, integer &innVp, double &inComTime,
		double *intimeSeries, double *infunSeries, double *ingradSeries, double *indistSeries, integer &inlengthSeries)
	{
		inx1 = x1; inf1 = f1; inngf0 = ngf0; inngf = ngf; initer = iter;
		innf = nf; inng = ng; innR = nR; innV = nV; innVp = nVp; inComTime = ComTime;
		intimeSeries = timeSeries; infunSeries = funSeries; ingradSeries = gradSeries; distSeries = indistSeries, inlengthSeries = lengthSeries;
	};

	void Solvers::PrintGenInfo(void)
	{
		printf("i:%d,f:%.3e,df/f:%.3e,", iter, f2, ((f1 - f2) / f2));

		if (nsubgf != -1)
			printf("nsubgf:%.3e,", nsubgf);

		printf("|gf|:%.3e,time:%.2e,", ngf, static_cast<double>(getTickCount() - starttime) / CLK_PS);

		if (subprobtimes != 0)
			printf("nsubprob:%d,", subprobtimes);

		printf("nf:%d,ng:%d,nR:%d,", nf, ng, nR);

		if (nV != 0)
			printf("nV(nVp):%d(%d),", nV, nVp);
	};

	void Solvers::PrintInfo(void)
	{
		printf("\n");
	};

	bool Solvers::IsStopped(void)
	{
		if (static_cast<double>(getTickCount() - starttime) / CLK_PS > TimeBound)
			return true;

		if (StopPtr != nullptr)
		{
			if (Prob->GetDomain()->GetIsIntrinsic())
			{
				const double * gf2space = gf2->GetSpace();
				if (gf2space != nullptr)
				{
					Vector *exgf2 = Prob->GetDomain()->GetEMPTYEXTR()->ConstructEmpty();
					Prob->GetDomain()->ObtainExtr(x2, gf2, exgf2);
					bool flag = StopPtr(x2, exgf2, f2, ngf, ngf0, Prob, this);
					delete exgf2;
					return flag;
				}
				return false;
			}
			else
			{
				return StopPtr(x2, gf2, f2, ngf, ngf0, Prob, this);
			}
		}

		if (Stop_Criterion == FUN_REL)
		{
			//if (fabs((f1 - f2) / f1) < Tolerance)
			//	numfref++;
			//else
			//	numfref = 0;
			//if (numfref > 4)
			//	return true;
			//return false;
			return ((fabs((f1 - f2) / (fabs(f1) + 1)) < Tolerance) && iter > 0);
		}
		else
		if (Stop_Criterion == GRAD_F)
			return ngf < Tolerance;
		else
		if (Stop_Criterion == GRAD_F_0)
			return (ngf / ngf0) < Tolerance;
		else
		if (Stop_Criterion == PSSUBGRAD)
		{
			if (gf2->GetSpace() == nullptr)
				return false;

			Variable *x1mx2 = x1->ConstructEmpty();
			Mani->VectorMinusVector(x1, x1, x2, x1mx2);
			if (sqrt(Mani->Metric(x1, x1mx2, x1mx2)) / (sqrt(Mani->Metric(x2, x2, x2)) + 1) >= Diffx)
			{
				gf2->CopyTo(gfs[0]);
				Currentlengthgfs = 1;
				idxgfs = 0;
			}
			else
			{
				if (Currentlengthgfs < Lengthgfs)
				{
					for (integer i = 0; i < Currentlengthgfs; i++)
					{
						Mani->VectorTransport(x1, eta2, x2, gfs[i], gfs[i]); nVp++;
					}
					gf2->CopyTo(gfs[Currentlengthgfs]);
					Currentlengthgfs++;
				}
				else
				{
					gf2->CopyTo(gfs[idxgfs]);
					idxgfs = (++idxgfs) % Lengthgfs;
					integer idx;
					for (integer i = idxgfs; i < idxgfs + Lengthgfs - 1; i++)
					{
						idx = i % Lengthgfs;
						Mani->VectorTransport(x1, eta2, x2, gfs[idx], gfs[idx]); nVp++;
					}
				}
			}
			delete x1mx2;

			if (Currentlengthgfs > 1)
				subprobtimes++;

			nsubgf = sqrt(MinPNormConHull(Mani, x2, gfs, Currentlengthgfs, nullptr, nullptr, 0));

			return fabs(nsubgf) < Tolerance;
		}

		printf("Error: Stopping Criterion is not specefic!\n");
		return true;
	};

	void Solvers::CheckParams(void)
	{
		std::string STOPCRITnames[STOPCRITLENGTH] = { "FUN_REL", "GRAD_F", "GRAD_F_0", "PSSUBGRAD" };
		std::string DEBUGnames[DEBUGLENGTH] = { "NOOUTPUT", "FINALRESULT", "ITERRESULT", "DETAILED" };
		char YES[] = "YES";
		char NO[] = "NO";
		char *status;
		printf("GENERAL PARAMETERS:\n");
		status = (Stop_Criterion >= 0 && Stop_Criterion < STOPCRITLENGTH) ? YES : NO;
		printf("Stop_Criterion:%15s[%s],\t", STOPCRITnames[Stop_Criterion].c_str(), status);
		status = (Tolerance > 0) ? YES : NO;
		printf("Tolerance     :%15g[%s]\n", Tolerance, status);
		status = (Max_Iteration > 0 && Max_Iteration >= Min_Iteration) ? YES : NO;
		printf("Max_Iteration :%15d[%s],\t", Max_Iteration, status);
		status = (Min_Iteration >= 0 && Min_Iteration <= Max_Iteration) ? YES : NO;
		printf("Min_Iteration :%15d[%s]\n", Min_Iteration, status);
		status = (OutputGap > 0) ? YES : NO;
		printf("OutputGap     :%15d[%s],\t", OutputGap, status);
		status = (Debug >= 0 && Debug < DEBUGLENGTH) ? YES : NO;
		printf("DEBUG         :%15s[%s]\n", DEBUGnames[Debug].c_str(), status);
		status = (Diffx > 0) ? YES : NO;
		printf("Diffx         :%15g[%s],\t", Diffx, status);
		status = (NumExtraGF > 0) ? YES : NO;
		printf("NumExtraGF    :%15d[%s]\n", NumExtraGF, status);
	};

	void Solvers::Run(void)
	{
		/*For partly smooth functions*/
		Lengthgfs = Mani->GetIntrDim() + NumExtraGF;
#ifdef MATLAB_MEX_FILE
		mxArray *lhs[1], *rhs[1];
		rhs[0] = mxCreateString("tic");
		mexCallMATLAB(0, lhs, 1, rhs, "feval");
#endif
		if (Stop_Criterion == PSSUBGRAD || SolverName == static_cast<std::string> ("LRBFGSLPSub") 
			|| SolverName == static_cast<std::string> ("RBFGSLPSub") || SolverName == static_cast<std::string> ("RGS"))
		{
			DeleteVectors(gfs, Lengthgfs);
			NewVectors(gfs, Lengthgfs);
		}
		if (SolverName == static_cast<std::string> ("RGS"))
		{
			DeleteVariables(Xs, Lengthgfs);
			NewVariables(Xs, Lengthgfs);
		}

		starttime = getTickCount();
		if (Debug >= ITERRESULT)
		{
			if (timeSeries != nullptr)
				delete[] timeSeries;
			timeSeries = new double[1 + Max_Iteration];
			if (funSeries == nullptr)
				delete[] funSeries;
			funSeries = new double[1 + Max_Iteration];
			if (gradSeries == nullptr)
				delete[] gradSeries;
			gradSeries = new double[1 + Max_Iteration];
			if (distSeries == nullptr)
				delete[] distSeries;
			distSeries = new double[1 + Max_Iteration];
		}
		if (Debug >= FINALRESULT)
			printf("=========================%s=========================\n", SolverName.c_str());
	};

	void Solvers::Initialization(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		SetProbX(prob, initialx, insoln);
		SetDefaultParams();
	};

	void Solvers::SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln)
	{
		const Vector *EMPTYETA;
		if (prob->GetDomain()->GetIsIntrinsic())
			EMPTYETA = prob->GetDomain()->GetEMPTYINTR();
		else
			EMPTYETA = prob->GetDomain()->GetEMPTYEXTR();

		Mani = prob->GetDomain();
		Prob = prob;
		if (insoln != nullptr)
		{
			soln = insoln->ConstructEmpty();
			insoln->CopyTo(soln);
		}
		else
		{
			soln = nullptr;
		}
		x1 = initialx->ConstructEmpty();

		initialx->CopyTo(x1);
		x2 = initialx->ConstructEmpty();
		gf1 = EMPTYETA->ConstructEmpty();
		gf2 = EMPTYETA->ConstructEmpty();

		eta1 = EMPTYETA->ConstructEmpty();
		eta2 = EMPTYETA->ConstructEmpty();
		zeta = EMPTYETA->ConstructEmpty();
	};

	void Solvers::SetDefaultParams()
	{
		nf = 0; ng = 0; nV = 0; nVp = 0; nR = 0; nH = 0; lengthSeries = 0;
		timeSeries = nullptr; funSeries = nullptr; gradSeries = nullptr; distSeries = nullptr;
		StopPtr = nullptr;

		Stop_Criterion = GRAD_F_0;
		TimeBound = 60 * 60 * 24 * 365;//one year;
		Tolerance = 1e-6;
		Max_Iteration = 500;
		Min_Iteration = 0;
		OutputGap = 1;
		Debug = ITERRESULT;
		nsubgf = -1;
		gfs = nullptr;
		NumExtraGF = 3;
		Currentlengthgfs = 0;
		idxgfs = 0;
		Diffx = 1e-6;
		subprobtimes = 0;
	};

	Solvers::~Solvers(void)
	{
		delete eta1;
		delete eta2;
		delete zeta;
		delete x1;
		delete x2;
		delete gf1;
		delete gf2;
		delete soln;
		if (Debug >= ITERRESULT)
		{
			if (timeSeries != nullptr)
				delete[] timeSeries;
			if (funSeries != nullptr)
				delete[] funSeries;
			if (gradSeries != nullptr)
				delete[] gradSeries;
			if (distSeries != nullptr)
				delete[] distSeries;
		}
		DeleteVectors(gfs, Lengthgfs);
	};

	void Solvers::NewVectors(Vector ** &Vs, integer l)
	{
		Vs = new Vector *[l];
		for (integer i = 0; i < l; i++)
			Vs[i] = gf1->ConstructEmpty();
	};

	void Solvers::DeleteVectors(Vector ** &Vs, integer l)
	{
		if (Vs != nullptr)
		{
			for (integer i = 0; i < l; i++)
				delete Vs[i];
			delete[] Vs;
		}
	};

	void Solvers::NewVariables(Vector ** &Xs, integer l)
	{
		Xs = new Vector *[l];
		for (integer i = 0; i < l; i++)
			Xs[i] = x1->ConstructEmpty();
	};

	void Solvers::DeleteVariables(Vector ** &Xs, integer l)
	{
		if (Xs != nullptr)
		{
			for (integer i = 0; i < l; i++)
				delete Xs[i];
			delete[] Xs;
		}
	};

	void Solvers::SetParams(PARAMSMAP params)
	{
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("Stop_Criterion"))
			{
				Stop_Criterion = static_cast<StopCrit> (static_cast<integer> (iter->second));
			}
			else
			if (iter->first == static_cast<std::string> ("Tolerance"))
			{
				Tolerance = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("TimeBound"))
			{
				TimeBound = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Max_Iteration"))
			{
				Max_Iteration = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("Min_Iteration"))
			{
				Min_Iteration = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("OutputGap"))
			{
				OutputGap = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("DEBUG"))
			{
				Debug = static_cast<DEBUGINFO> (static_cast<integer> (iter->second));
			}
			else
			if (iter->first == static_cast<std::string> ("Diffx"))
			{
				Diffx = static_cast<double> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("NumExtraGF"))
			{
				NumExtraGF = static_cast<integer> (iter->second);
			}
		}
	};
}; /*end of ROPTLIB namespace*/
