
#include "Solvers/QuasiNewton.h"

/*Define the namespace*/
namespace ROPTLIB{

	void QuasiNewton::HvRBroydenFamily(Vector *v, Vector *result)
	{
		Mani->LinearOPEEta(x1, H, v, result); nH++;
	};

	void QuasiNewton::UpdateDataRBroydenFamily(void)
	{
		double yHy;
		Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
		Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
		betay = Mani->Beta(x1, eta2);
		//Mani->VectorLinearCombination(x2, 1.0 / betay, gf2, -1.0, zeta, y);
		Mani->scalarVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);
		inpsy = Mani->Metric(x2, s, y);
		if (isconvex && iter == 1 && inpsy > 0)
			H->ScaledIdOPE(inpsy / Mani->Metric(x2, y, y));
		Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
		inpss = Mani->Metric(x2, s, s);
		if (inpsy / inpss >= nu * pow(ngf, mu) && (ngf / ngf0 < 1e-3 ||
			(inpss > std::numeric_limits<double>::epsilon() && inpsy > std::numeric_limits<double>::epsilon())))
		{
			Mani->LinearOPEEta(x2, tildeH, y, zeta); // zeta = tildeH y
			yHy = Mani->Metric(x2, y, zeta);
			Mani->VectorLinearCombination(x2, 1.0 / inpsy, s, -1.0 / yHy, zeta, u);
			phic = Phi(x2, y, s, tildeH, inpsy, yHy, u);
			Mani->HaddScaledRank1OPE(x2, tildeH, -1.0 / yHy, zeta, zeta, H);
			Mani->HaddScaledRank1OPE(x2, H, 1.0 / inpsy, s, s, H);
			Mani->HaddScaledRank1OPE(x2, H, phic * yHy, u, u, H);
			isupdated = true;
		}
		else
		{
			isupdated = false;
			tildeH->CopyTo(H);
		}
	};

	double QuasiNewton::Phi(Variable *x2, Vector *y, Vector *s, LinearOPE *tildeH, double inpsy, double yHy, Vector *u)
	{
		return 1;
	};

	void QuasiNewton::HvRBFGSSub(Vector *v, Vector *result)
	{
		//v->CopyTo(result);
		//return;//---
		Mani->LinearOPEEta(x1, H, v, result); nH++;
	};

	void QuasiNewton::UpdateDataRBFGSSub(void)
	{
		Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
		Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
		betay = Mani->Beta(x1, eta2);
		Mani->scalarVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);
		inpyy = Mani->Metric(x2, y, y);
		inpsy = Mani->Metric(x2, s, y);
		double tmp = 1.0 / lambdaUpper - inpsy / inpyy;
		tmp = (tmp > 0) ? tmp : 0;
		Mani->scalarVectorAddVector(x2, tmp, y, s, s);
		inpsy = Mani->Metric(x2, s, y);

		if (isconvex && iter == 1 && inpsy > 0)
			H->ScaledIdOPE(inpsy / inpyy);

		Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
		inpss = Mani->Metric(x2, s, s);
		if (inpsy / inpss > lambdaLower)
		{
			Mani->LinearOPEEta(x2, tildeH, y, zeta); // zeta = tildeH y
			Mani->HaddScaledRank1OPE(x2, tildeH, -1.0 / inpsy, s, zeta, H);
			Mani->LinearOPEEta(x2, H, y, zeta); // zeta = H y

			Mani->HaddScaledRank1OPE(x2, H, -1.0 / inpsy, zeta, s, H);
			Mani->HaddScaledRank1OPE(x2, H, 1.0 / inpsy, s, s, H);
			isupdated = true;
		}
		else
		{
			isupdated = false;
			H->ScaledIdOPE(1);
		}
	};

	void QuasiNewton::HvRBFGS(Vector *v, Vector *result)
	{
		Mani->LinearOPEEta(x1, H, v, result); nH++;
	};

	void QuasiNewton::UpdateDataRBFGS(void)
	{
		Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
		Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
		betay = Mani->Beta(x1, eta2);
		Mani->scalarVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);
		inpsy = Mani->Metric(x2, s, y);
		if (isconvex && iter == 1 && inpsy > 0)
			H->ScaledIdOPE(inpsy / Mani->Metric(x2, y, y));

		Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
		inpss = Mani->Metric(x2, s, s);
		if ((inpsy / inpss >= nu * pow(ngf, mu) && (ngf / ngf0 < 1e-3 ||
			(inpss > std::numeric_limits<double>::epsilon() && inpsy > std::numeric_limits<double>::epsilon()))) 
			|| (Stop_Criterion == PSSUBGRAD))
		{
			Mani->LinearOPEEta(x2, tildeH, y, zeta); // zeta = tildeH y
			Mani->HaddScaledRank1OPE(x2, tildeH, -1.0 / inpsy, s, zeta, H);
			Mani->LinearOPEEta(x2, H, y, zeta); // zeta = H y

			Mani->HaddScaledRank1OPE(x2, H, -1.0 / inpsy, zeta, s, H);
			Mani->HaddScaledRank1OPE(x2, H, 1.0 / inpsy, s, s, H);
			isupdated = true;
		}
		else
		{
			isupdated = false;
			tildeH->CopyTo(H);
#if defined(TESTELASTICCURVESRO) //-- || defined(TESTORTHBOUNDINGBOX)
			H->ScaledIdOPE(1);
#endif
		}
	};


	void QuasiNewton::HvLRBFGSSub(Vector *v, Vector *result)
	{
		double *xi = new double[Currentlength];
		double omega;
		integer idx;
		v->CopyTo(result);
		for (integer i = Currentlength - 1; i >= 0; i--)
		{
			idx = (beginidx + i) % LengthSY;
			xi[idx] = RHO[idx] * Mani->Metric(x1, S[idx], result);
			Mani->scalarVectorAddVector(x1, -xi[idx], Y[idx], result, result);
		}
		Mani->ScaleTimesVector(x1, gamma, result, result);
		for (integer i = 0; i < Currentlength; i++)
		{
			idx = (beginidx + i) % LengthSY;
			omega = RHO[idx] * Mani->Metric(x1, Y[idx], result);
			Mani->scalarVectorAddVector(x1, xi[idx] - omega, S[idx], result, result);
		}

		delete[] xi;
	};

	void QuasiNewton::UpdateDataLRBFGSSub(void)
	{
		Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
		Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
		betay = Mani->Beta(x1, eta2);
		Mani->scalarVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);

		inpyy = Mani->Metric(x2, y, y);
		inpsy = Mani->Metric(x2, s, y);
		double tmp = 1.0 / lambdaUpper - inpsy / inpyy;
		tmp = (tmp > 0) ? tmp : 0;
		Mani->scalarVectorAddVector(x2, tmp, y, s, s);
		inpsy = Mani->Metric(x2, s, y);
		inpss = Mani->Metric(x2, s, s);
		rho = 1.0 / inpsy;
		if (inpsy / inpss >= lambdaLower)
		{
			gamma = 1;
			//gamma = inpsy / inpyy; /*Suggested in NW2006*/
			//gamma = inpss / inpsy; /*BB stepsize*/
			if (Currentlength < LengthSY)
			{
				y->CopyTo(Y[Currentlength]);
				s->CopyTo(S[Currentlength]);
				RHO[Currentlength] = rho;
				for (integer i = 0; i < Currentlength; i++)
				{
					Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]); nVp++;
					Mani->VectorTransport(x1, eta2, x2, S[i], S[i]); nVp++;
				}
				Currentlength++;
			}
			else
				if (LengthSY > 0)
				{
					integer idx;
					y->CopyTo(Y[beginidx]);
					s->CopyTo(S[beginidx]);
					RHO[beginidx] = rho;
					beginidx = (++beginidx) % LengthSY;
					for (integer i = beginidx; i < beginidx + LengthSY - 1; i++)
					{
						idx = i % LengthSY;
						Mani->VectorTransport(x1, eta2, x2, Y[idx], Y[idx]); nVp++;
						Mani->VectorTransport(x1, eta2, x2, S[idx], S[idx]); nVp++;
					}
				}
			isupdated = true;
		}
		else
		{
			Currentlength = 0;
			beginidx = 0;
			//for (integer i = 0; i < Currentlength; i++)
			//{
			//	Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]); nVp++;
			//	Mani->VectorTransport(x1, eta2, x2, S[i], S[i]); nVp++;
			//}
			isupdated = false;
		}
	};

	void QuasiNewton::HvLRBFGS(Vector *v, Vector *result)
	{
		double *xi = new double[Currentlength];
		double omega;
		integer idx;
		v->CopyTo(result);
		for (integer i = Currentlength - 1; i >= 0; i--)
		{
			idx = (beginidx + i) % LengthSY;
			xi[idx] = RHO[idx] * Mani->Metric(x1, S[idx], result);
			Mani->scalarVectorAddVector(x1, -xi[idx], Y[idx], result, result);
		}
		Mani->ScaleTimesVector(x1, gamma, result, result);
		for (integer i = 0; i < Currentlength; i++)
		{
			idx = (beginidx + i) % LengthSY;
			omega = RHO[idx] * Mani->Metric(x1, Y[idx], result);
			Mani->scalarVectorAddVector(x1, xi[idx] - omega, S[idx], result, result);
		}
		delete[] xi;
	};

	double QuasiNewton::InitialHessian(double inpss, double inpsy, double inpyy)
	{ /*Suggested in NW2006*/
		return inpsy / inpyy;
	};

	void QuasiNewton::UpdateDataLRBFGS(void)
	{
#ifdef TESTEUCPOSSPCD
		/*For sparse coding problem, the domain is not a manifold.
		If the active changes, then restart the LRBFGS.*/
		const double *x1ptr = x1->ObtainReadData();
		const double *x2ptr = x2->ObtainReadData();
		for (integer i = 0; i < x1->Getlength(); i++)
		{
			if (x1ptr[i] != 0 && x2ptr[i] <= Tolerance)
			{
				Currentlength = 0;
				beginidx = 0;
				if(Debug >= ITERRESULT && iter % OutputGap == 0)
				{
					printf("%d: Restart LRBFGS.\n", iter);
				}
				return;
			}
		}
#endif
		Mani->VectorTransport(x1, eta2, x2, eta2, s); nV++;
		Mani->VectorTransport(x1, eta2, x2, gf1, zeta); nVp++;
		betay = Mani->Beta(x1, eta2);
		Mani->scalarVectorMinusVector(x2, 1.0 / betay, gf2, zeta, y);
		inpsy = Mani->Metric(x2, s, y);
		inpss = Mani->Metric(x2, s, s);
		inpyy = Mani->Metric(x2, y, y);
		rho = 1.0 / inpsy;
		if (inpsy / inpss >= nu * pow(ngf, mu) && (ngf / ngf0 < 1e-3 ||
			(inpss > std::numeric_limits<double>::epsilon() && inpsy > std::numeric_limits<double>::epsilon())))
		{
			gamma = InitialHessian(inpss, inpsy, inpyy);
			if (Currentlength < LengthSY)
			{
				y->CopyTo(Y[Currentlength]);
				s->CopyTo(S[Currentlength]);
				RHO[Currentlength] = rho;
				for (integer i = 0; i < Currentlength; i++)
				{
					Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]); nVp++;
					Mani->VectorTransport(x1, eta2, x2, S[i], S[i]); nVp++;
				}
				Currentlength++;
			}
			else
				if (LengthSY > 0)
				{
					integer idx;
					y->CopyTo(Y[beginidx]);
					s->CopyTo(S[beginidx]);
					RHO[beginidx] = rho;
					beginidx = (++beginidx) % LengthSY;
					for (integer i = beginidx; i < beginidx + LengthSY - 1; i++)
					{
						idx = i % LengthSY;
						Mani->VectorTransport(x1, eta2, x2, Y[idx], Y[idx]); nVp++;
						Mani->VectorTransport(x1, eta2, x2, S[idx], S[idx]); nVp++;
					}
				}
			isupdated = true;
		}
		else
		{
			for (integer i = 0; i < Currentlength; i++)
			{
				Mani->VectorTransport(x1, eta2, x2, Y[i], Y[i]); nVp++;
				Mani->VectorTransport(x1, eta2, x2, S[i], S[i]); nVp++;
			}
			isupdated = false;
		}
	};

	void QuasiNewton::HvRWRBFGS(Vector *v, Vector *result)
	{
		Mani->LinearOPEEta(x1, H, v, result); nH++;
	};

	void QuasiNewton::UpdateDataRWRBFGS(void)
	{
		eta2->CopyTo(s);
		Mani->coTangentVector(x1, eta2, x2, gf2, y); nV++;
		Mani->VectorMinusVector(x1, y, gf1, y);
		inpsy = Mani->Metric(x1, s, y);
		if (isconvex && iter == 1 && inpsy > 0)
			H->ScaledIdOPE(inpsy / Mani->Metric(x1, y, y));
		inpss = Mani->Metric(x1, s, s);
		if (inpsy / inpss >= nu * pow(ngf, mu) && inpss > std::numeric_limits<double>::epsilon()
			&& inpsy > std::numeric_limits<double>::epsilon())
		{
			Mani->LinearOPEEta(x2, H, y, zeta); // zeta = H y
			Mani->HaddScaledRank1OPE(x2, H, -1.0 / inpsy, s, zeta, H);
			Mani->LinearOPEEta(x2, H, y, zeta); // zeta = H y

			Mani->HaddScaledRank1OPE(x2, H, -1.0 / inpsy, zeta, s, H);
			Mani->HaddScaledRank1OPE(x2, H, 1.0 / inpsy, s, s, H);
			Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
			tildeH->CopyTo(H);
			isupdated = true;
		}
		else
		{
			isupdated = false;
			Mani->TranHInvTran(x1, eta2, x2, H, tildeH);
			tildeH->CopyTo(H);
		}
	};

	void QuasiNewton::HvRTRSR1(Vector *v, Vector *result)
	{
		Mani->LinearOPEEta(x1, B, v, result);
	};

	void QuasiNewton::UpdateDataRTRSR1(void)
	{
		double denorminator, norm2ymBs;
		double mintolsq = std::numeric_limits<double>::epsilon();
		Prob->Grad(x2, gf2); ng++;
		eta2->CopyTo(s);
		Mani->InverseVectorTransport(x1, eta2, x2, gf2, eta1); nV++;
		Mani->VectorMinusVector(x1, eta1, gf1, y);
		Mani->VectorMinusVector(x1, y, zeta, zeta);
		denorminator = Mani->Metric(x1, s, zeta);
		if (isconvex && iter == 1 && Mani->Metric(x1, s, y) > 0)
			B->ScaledIdOPE(Mani->Metric(x1, y, y) / Mani->Metric(x1, s, y));
		inpss = Mani->Metric(x1, s, s);
		norm2ymBs = Mani->Metric(x1, zeta, zeta);
		if (denorminator * denorminator >= mintolsq * inpss * norm2ymBs && (norm2ymBs >= mintolsq || ngf / ngf0 < 1e-3))
		{
			Mani->HaddScaledRank1OPE(x1, B, 1.0 / denorminator, zeta, zeta, B);
			isupdated = true;
		}
		else
		{
			isupdated = false;
		}
	};

	void QuasiNewton::HvLRTRSR1(Vector *Eta, Vector *result)
	{
		/* This function makes use of SS, SY and gamma to evaluate the action of Hessian approximation [HAG2014, (64)].
		[HAG2014]: W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trustregion method.
		Mathematical Programming, 150(2):179?16, February 2015.

		SS is the Q in (46), SY is the P in (46), PMGQ is the P - gamma Q in (46).
		*/
		integer idx;
		double *v = new double[Currentlength];

		/*if S and Y has been updated in function: UpdateData(void), then PMGQ is recomputed and
		the LU decomposition of PMGQ is also recomputed. The LU decomposition is used to evaluate
		the action of PMGQ^{-1}. */
		if (ischangedSandY)
		{
			for (integer i = 0; i < Currentlength; i++)
			{
				idx = (i + beginidx) % LengthSY;
				Mani->scalarVectorAddVector(x1, -gamma, S[idx], Y[idx], YMGS[i]);
			}
			for (integer i = 0; i < Currentlength; i++)
			{
				for (integer j = 0; j < Currentlength; j++)
				{
					PMGQ[i + j * Currentlength] = SY[i + j * LengthSY] - gamma * SS[i + j * LengthSY];
				}
			}
			if (Currentlength > 0)
			{
				// compute LU
				integer info, CurLen = Currentlength;
				// LU decomposion for PMGQ, PMGQ = P * L * U, L and U are stored in PMGQ, the permutation matrix is in P
				// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
				dgetrf_(&CurLen, &CurLen, PMGQ, &CurLen, P, &info);
				ischangedSandY = false;
			}
		}

		for (integer i = 0; i < Currentlength; i++)
			v[i] = Mani->Metric(x1, YMGS[i], Eta);

		if (Currentlength > 0)
		{
			char *trans = const_cast<char *> ("n");
			integer info, one = 1, CurLen = Currentlength;
			// solve linear system: PMGQ * X = v using the LU decomposition results from dgetrf, then solution is stored in v.
			// details: http://www.netlib.org/lapack/explore-html/d6/d49/dgetrs_8f.html
			dgetrs_(trans, &CurLen, &one, PMGQ, &CurLen, P, v, &CurLen, &info);
		}

		Mani->ScaleTimesVector(x1, gamma, Eta, result);
		for (integer i = 0; i < Currentlength; i++)
		{
			Mani->scalarVectorAddVector(x1, v[i], YMGS[i], result, result);
		}

		delete[] v;
	};

	void QuasiNewton::UpdateDataLRTRSR1(void)
	{
		double denorminator, norm2ymBs;
		double mintolsq = std::numeric_limits<double>::epsilon();
		double mintol = sqrt(mintolsq);
		Prob->Grad(x2, gf2); ng++;
		eta2->CopyTo(s);
		Mani->InverseVectorTransport(x1, eta2, x2, gf2, eta1); nV++;
		Mani->VectorMinusVector(x1, eta1, gf1, y);
		Mani->VectorMinusVector(x1, y, zeta, zeta);
		denorminator = Mani->Metric(x1, s, zeta);
		inpss = Mani->Metric(x1, s, s);
		norm2ymBs = Mani->Metric(x1, zeta, zeta);
		if (iter == 0) // This is for the robustness when the cost function is quadratic 
		{			   // and its Hessian is identity everywhere.
			inpsy = Mani->Metric(x1, s, y);
			inpyy = Mani->Metric(x1, y, y);
			gamma = inpyy / inpsy;
		}
		if (denorminator * denorminator >= mintolsq * inpss * norm2ymBs && (norm2ymBs >= mintolsq || ngf / ngf0 < 1e-3)
			&& (iter != 0 || fabs(gamma - inpsy / inpss) > mintol)) // This is for the robustness when the cost 
			// function is quadratic and its Hessian is identity everywhere.
		{
			inpsy = Mani->Metric(x1, s, y);
			inpyy = Mani->Metric(x1, y, y);
			gamma = inpyy / inpsy;

			/*if s and y are accepted, then S and Y need to be updated. It follows that the matrices SY and SS need to be update.*/
			if (Currentlength < LengthSY)
			{
				s->CopyTo(S[Currentlength]);
				y->CopyTo(Y[Currentlength]);
				SS[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], S[Currentlength]);
				SY[Currentlength + Currentlength * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[Currentlength]);
				for (integer i = 0; i < Currentlength; i++)
				{
					SS[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], S[i]);
					SS[i + Currentlength * LengthSY] = SS[Currentlength + i * LengthSY];
					SY[Currentlength + i * LengthSY] = Mani->Metric(x1, S[Currentlength], Y[i]);
					SY[i + Currentlength * LengthSY] = SY[Currentlength + i * LengthSY];
				}
				Currentlength++;
			}
			else
			{
				s->CopyTo(S[beginidx]);
				y->CopyTo(Y[beginidx]);
				for (integer i = 0; i < LengthSY - 1; i++)
				{
					for (integer j = 0; j < LengthSY - 1; j++)
					{
						SS[i + j * LengthSY] = SS[i + 1 + (j + 1) * LengthSY];
						SY[i + j * LengthSY] = SY[i + 1 + (j + 1) * LengthSY];
					}
				}
				SS[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], S[beginidx]);
				SY[LengthSY * LengthSY - 1] = Mani->Metric(x1, S[beginidx], Y[beginidx]);
				integer idx = 0;
				for (integer i = 0; i < LengthSY - 1; i++)
				{
					idx = (i + beginidx + 1) % LengthSY;
					SS[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, S[idx], S[beginidx]);
					SS[LengthSY - 1 + i * LengthSY] = SS[i + (LengthSY - 1) * LengthSY];
					SY[i + (LengthSY - 1) * LengthSY] = Mani->Metric(x1, Y[idx], S[beginidx]);
					SY[LengthSY - 1 + i * LengthSY] = SY[i + (LengthSY - 1) * LengthSY];
				}
				beginidx = (++beginidx) % LengthSY;
			}
			isupdated = true;
			ischangedSandY = true;
		}
		else
		{
			isupdated = false;
		}
	};

	void QuasiNewton::SetParams(PARAMSMAP params)
	{
		Solvers::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("isconvex"))
			{
				isconvex = ((static_cast<integer> (iter->second)) != 0);
			}
			else
			if (iter->first == static_cast<std::string> ("LengthSY"))
			{
				LengthSY = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("nu"))
			{
				nu = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("mu"))
			{
				mu = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("lambdaLower"))
			{
				lambdaLower = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("lambdaUpper"))
			{
				lambdaUpper = iter->second;
			}
			else
			if (iter->first == static_cast<std::string> ("Num_pre_BB"))
			{
				Num_pre_BB = static_cast<integer> (iter->second);
			}
			else
			if (iter->first == static_cast<std::string> ("BBratio"))
			{
				BBratio = static_cast<double> (iter->second);
			}
		}
	};
}; /*end of ROPTLIB namespace*/
