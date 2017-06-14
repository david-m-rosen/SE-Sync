
#include "Manifolds/Sphere/Sphere.h"

/*Define the namespace*/
namespace ROPTLIB{

	Sphere::Sphere(integer inn) :Stiefel(inn, 1)
	{
		name.assign("Sphere");
		delete EMPTYEXTR;
		delete EMPTYINTR;
		EMPTYEXTR = new SphereVector(n);
		EMPTYINTR = new SphereVector(IntrinsicDim);
	};

	Sphere::~Sphere(void)
	{
	};

	// choose qf retraction, parallelization and intrinsic approach and no householder reflections
	void Sphere::ChooseSphereParamsSet1(void)
	{
		metric = EUCLIDEAN;
		retraction = QF;
		VecTran = PARALLELIZATION;
		IsIntrApproach = true;
		HasHHR = false;
		UpdBetaAlone = false;
		HasLockCon = false;
	};
	// choose exponential map, parallel translation and extrinsic approach and no householder reflections
	// Even though the Householder reflections are not used, the locking condition is satisfied.
	void Sphere::ChooseSphereParamsSet2(void)
	{
		metric = EUCLIDEAN;
		retraction = EXP;
		VecTran = PARALLELTRANSLATION;
		IsIntrApproach = false;
		HasHHR = false;
		UpdBetaAlone = false;
		HasLockCon = true;
	};

	// choose qf, parallel translation and extrinsic approach and no householder reflections
	// The locking conidition is not satisfied
	void Sphere::ChooseSphereParamsSet3(void)
	{
		metric = EUCLIDEAN;
		retraction = QF;
		VecTran = PARALLELTRANSLATION;
		IsIntrApproach = false;
		HasHHR = false;
		UpdBetaAlone = false;
		HasLockCon = false;
	};

	// choose qf, parallel translation and extrinsic approach and no householder reflections
	// Beta \neq 1 is used and the locking conidition is satisfied
	void Sphere::ChooseSphereParamsSet4(void)
	{
		metric = EUCLIDEAN;
		retraction = QF;
		VecTran = PARALLELTRANSLATION;
		IsIntrApproach = false;
		HasHHR = false;
		UpdBetaAlone = true;
		HasLockCon = true;
	};

	/* choose proximal mapping, parallelization and intrinsic approach and no householder reflections
	the locking conidition is not satisfied*/
	void Sphere::ChooseSphereParamsSet5(void)
	{
		metric = EUCLIDEAN;
		retraction = PROXSTIE;
		VecTran = PARALLELIZATION;
		IsIntrApproach = true;
		HasHHR = false;
		UpdBetaAlone = false;
		HasLockCon = false;
	};

	void Sphere::ProxRetraction(Variable *x, Vector *etax, Variable *result, double instepsize) const
	{
		if (!GetIsIntrinsic())
		{
			VectorAddVector(x, x, etax, result);
		}
		else
		{
			Vector *exetax = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(x, etax, exetax);
			VectorAddVector(x, x, exetax, result);
			delete exetax;
		}
		
		double *resultptr = result->ObtainWritePartialData();
		double maxv = 0;
		integer maxi = 0;
		for (integer i = 0; i < result->Getlength(); i++)
		{
			if (fabs(resultptr[i]) > fabs(maxv))
			{
				maxv = resultptr[i];
				maxi = i;
			}
		}
		for (integer i = 0; i < result->Getlength(); i++)
		{
			resultptr[i] = (fabs(resultptr[i]) < instepsize) ? 0 : ((resultptr[i] > instepsize) ? resultptr[i] - instepsize : resultptr[i] + instepsize);
		}
		double rnorm = sqrt(Metric(x, result, result));
		if (rnorm != 0)
			ScaleTimesVector(x, 1.0 / rnorm, result, result);
		else
		{
			for (integer i = 0; i < result->Getlength(); i++)
			{
				resultptr[i] = 0;
			}
			resultptr[maxi] = (maxv >= 0) ? 1 : -1;
		}
	};

	void Sphere::ExpRetraction(Variable *x, Vector *etax, Variable *result) const
	{
		double normetax = sqrt(Metric(x, etax, etax));
		VectorLinearCombination(x, cos(normetax), x, sin(normetax) / normetax, etax, result);
		double normresult = sqrt(Metric(x, result, result));
		ScaleTimesVector(x, 1.0 / normresult, result, result);
	};

	void Sphere::ExpcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		double xiytx = Metric(x, x, xiy);
		double xiytetax = Metric(x, xiy, etax);
		double normetax = sqrt(Metric(x, etax, etax));
		double sinnormetax = sin(normetax);
		double cosnormetax = cos(normetax);
		VectorLinearCombination(x, sinnormetax / normetax, xiy, (xiytetax * cosnormetax / normetax
			- xiytx * sinnormetax - xiytetax * sinnormetax / normetax / normetax) / normetax, etax, result);
		Projection(x, result, result);
	};

	void Sphere::ExpDiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		double etaxtxix = Metric(x, etax, xix);
		double normetax = sqrt(Metric(x, etax, etax));
		double sinnormetax = sin(normetax);
		VectorLinearCombination(x, -sinnormetax * etaxtxix / normetax, x, sinnormetax / normetax, xix, result);
		scalarVectorAddVector(x, (cos(normetax) - sinnormetax / normetax) * etaxtxix / normetax / normetax, etax, result, result);
	};

	void Sphere::ExpVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		Vector *xpy = x->ConstructEmpty();
		VectorAddVector(x, x, y, xpy);
		double xynormsq = Metric(x, xpy, xpy);
		double xixty = Metric(x, xix, y);
		scalarVectorAddVector(x, -2.0 * xixty / xynormsq, xpy, xix, result);
		delete xpy;
	};

	void Sphere::ExpInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		Vector *xpy = x->ConstructEmpty();
		VectorAddVector(x, x, y, xpy);
		double xynormsq = Metric(x, xpy, xpy);
		double xiytx = Metric(x, xiy, x);
		scalarVectorAddVector(x, -2.0 * xiytx / xynormsq, xpy, xiy, result);
		delete xpy;
	};

	void Sphere::ExpHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		Vector *xpy = etax->ConstructEmpty();
		VectorAddVector(x, x, y, xpy);
		integer ell = Hx->Getsize()[0];
		integer length = etax->Getlength();
		const double *M = Hx->ObtainReadData();
		double *Hxpy = new double[ell];
		const double *xpyTV = xpy->ObtainReadData();

		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		integer inc = 1, N = ell;
		// Hxpy <- M(: start : start + length - 1) * xpyTV, details: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
		dgemv_(transn, &N, &length, &one, const_cast<double *> (M + start * N), &N, const_cast<double *> (xpyTV), &inc, &zero, Hxpy, &inc);

		double scalar = -2.0 / Metric(x, xpy, xpy);
		Hx->CopyTo(result);
		const double *xv = x->ObtainReadData();
		double *resultL = result->ObtainWritePartialData();
		// resultL(:, start : start + length - 1) <- scalar * Hxpy * xv^T + resultL(:, start : start + length - 1),
		// details: http://www.netlib.org/lapack/explore-html/dc/da8/dger_8f.html
		dger_(&length, &N, &scalar, Hxpy, &inc, const_cast<double *> (xv), &inc, resultL + start * N, &N);

		delete xpy;
		delete[] Hxpy;
	};

	void Sphere::ExpTranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		Vector *xpy = etax->ConstructEmpty();
		VectorAddVector(x, x, y, xpy);
		integer ell = Hx->Getsize()[0];
		integer length = etax->Getlength();
		const double *M = Hx->ObtainReadData();
		double *Hty = new double[ell];
		const double *yv = y->ObtainReadData();

		char *transt = const_cast<char *> ("t");
		double one = 1, zero = 0;
		integer inc = 1, N = ell;
		// Hty <- M(start : start + length - 1, :)^T * yv, details: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
		dgemv_(transt, &length, &N, &one, const_cast<double *> (M + start), &N, const_cast<double *> (yv), &inc, &zero, Hty, &inc);

		double scalar = -2.0 / Metric(x, xpy, xpy);
		const double *xpyTV = xpy->ObtainReadData();
		Hx->CopyTo(result);
		double *resultL = result->ObtainWritePartialData();
		// resultL(start : start + length - 1, :) <- scalar * xpyTV * Hty^T + resultL(start : start + length - 1, :),
		// details: http://www.netlib.org/lapack/explore-html/dc/da8/dger_8f.html
		dger_(&length, &N, &scalar, const_cast<double *> (xpyTV), &inc, Hty, &inc, resultL + start, &N);

		delete xpy;
		delete[] Hty;
	};

	void Sphere::ExpTranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
	{
		ExpHInvTran(x, etax, y, Hx, 0, etax->Getlength(), result);
		ExpTranH(x, etax, y, result, 0, etax->Getlength(), result);
	};

	void Sphere::Retraction(Variable *x, Vector *etax, Variable *result) const
	{
		Retraction(x, etax, result, 1.0);
	};

	void Sphere::Retraction(Variable *x, Vector *etax, Variable *result, double instepsize) const
	{
		if (retraction == EXP)
		{
			ExpRetraction(x, etax, result);
			return;
		}
		if (retraction == PROXSTIE)
		{
			ProxRetraction(x, etax, result, instepsize);
			return;
		}

		Stiefel::Retraction(x, etax, result);
	};

	void Sphere::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (retraction == EXP)
		{
			ExpcoTangentVector(x, etax, y, xiy, result);
			return;
		}
		Stiefel::coTangentVector(x, etax, y, xiy, result);
	};

	void Sphere::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (retraction == EXP)
		{
			ExpDiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
			return;
		}
		Stiefel::DiffRetraction(x, etax, y, xix, result, IsEtaXiSameDir);
	};

	void Sphere::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		if (VecTran == PARALLELTRANSLATION)
		{
			ExpVectorTransport(x, etax, y, xix, result);
			return;
		}
		Stiefel::VectorTransport(x, etax, y, xix, result);
	};

	void Sphere::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (VecTran == PARALLELTRANSLATION)
		{
			ExpInverseVectorTransport(x, etax, y, xiy, result);
			return;
		}
		Stiefel::InverseVectorTransport(x, etax, y, xiy, result);
	};

	void Sphere::HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == PARALLELTRANSLATION)
		{
			ExpHInvTran(x, etax, y, Hx, start, end, result);
			return;
		}
		Stiefel::HInvTran(x, etax, y, Hx, start, end, result);
	};

	void Sphere::TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == PARALLELTRANSLATION)
		{
			ExpTranH(x, etax, y, Hx, start, end, result);
			return;
		}
		Stiefel::TranH(x, etax, y, Hx, start, end, result);
	};

	void Sphere::TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
	{
		if (VecTran == PARALLELTRANSLATION)
		{
			ExpTranHInvTran(x, etax, y, Hx, result);
			return;
		}
		Stiefel::TranHInvTran(x, etax, y, Hx, result);
	};

	void Sphere::SetParams(PARAMSMAP params)
	{
		Stiefel::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("ParamSet"))
			{
				switch (static_cast<integer> (iter->second))
				{
				case 1:
					ChooseSphereParamsSet1();
					break;
				case 2:
					ChooseSphereParamsSet2();
					break;
				case 3:
					ChooseSphereParamsSet3();
					break;
				case 4:
					ChooseSphereParamsSet4();
					break;
				case 5:
					ChooseSphereParamsSet5();
					break;
				default:
					break;
				}
			}
		}
	};
}; /*end of ROPTLIB namespace*/
