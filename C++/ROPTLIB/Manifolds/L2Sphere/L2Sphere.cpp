
#include "Manifolds/L2Sphere/L2Sphere.h"

/*Define the namespace*/
namespace ROPTLIB{

	L2Sphere::L2Sphere(integer inn)
	{
		// public parameter
		HasHHR = false;

		// Only some combinations exist.
		metric = TRAPEZOID;
		retraction = NORMALIZED;
		VecTran = L2SPARALLELTRANSLATION;
		IsIntrApproach = false;
		UpdBetaAlone = false;

		// Status of locking condition
		HasLockCon = false;

		// Fixed parameters
		n = inn;
		ExtrinsicDim = n;
		IntrinsicDim = n - 1;
		name.assign("L2Sphere");

		EMPTYEXTR = new L2SphereVector(n);
		EMPTYINTR = new L2SphereVector(IntrinsicDim);
	};

	L2Sphere::~L2Sphere()
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	double L2Sphere::Metric(Variable *x, Vector *etax, Vector *xix) const
	{ //Trapezoidal rule
		const double *etaxTV = etax->ObtainReadData();
		const double *xixTV = xix->ObtainReadData();
		integer inc = 1;
		// In Matlab environment, the following library produces different results from different runs, which
		// I don't know why. However, it can be fixed by the next 3 lines of codes, which is slower.

		// output etaxTV^T xixTV, details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		double result = ddot_(&n, const_cast<double *> (etaxTV), &inc, const_cast<double *> (xixTV), &inc);
		//double result = 0;
		//for (integer i = 0; i < n; i++)
		//	result += etaxTV[i] * xixTV[i];
		result -= etaxTV[0] * xixTV[0] / 2;
		result -= etaxTV[n - 1] * xixTV[n - 1] / 2;
		return result / (n - 1);
	};

	void L2Sphere::Projection(Variable *x, Vector *v, Vector *result) const
	{
		const double *xl = x->ObtainReadData();
		double nume = Metric(x, x, v);
		scalarVectorAddVector(x, -nume, x, v, result);
	};

	void L2Sphere::Retraction(Variable *x, Vector *etax, Variable *result) const
	{// exponential mapping
		double norm = sqrt(Metric(x, etax, etax));
		if (norm < std::numeric_limits<double>::epsilon())
			ScaleTimesVector(x, cos(norm), x, result);
		else
			VectorLinearCombination(x, cos(norm), x, sin(norm) / norm, etax, result);

		norm = sqrt(this->Metric(x, result, result));
		this->ScaleTimesVector(x, 1.0 / norm, result, result);
	};

	void L2Sphere::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		xiy->CopyTo(result);
		printf("The cotangent vector has not been implemented!\n");
	};

	void L2Sphere::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsEtaXiSameDir)
		{
			VectorTransport(x, etax, y, xix, result);

			if (IsEtaXiSameDir && (HasHHR || UpdBetaAlone))
			{
				const double *etaxTV = etax->ObtainReadData();
				const double *xixTV = xix->ObtainReadData();
				double EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
				SharedSpace *beta = new SharedSpace(1, 3);
				double *betav = beta->ObtainWriteEntireData();
				betav[0] = sqrt(Metric(x, etax, etax) / Metric(x, result, result)) / EtatoXi;
				betav[1] = Metric(x, etax, etax);
				betav[2] = Metric(x, result, result) * EtatoXi * EtatoXi;
				etax->AddToTempData("beta", beta);

				if (HasHHR)
				{
					Vector *TReta = result->ConstructEmpty();
					result->CopyTo(TReta);
					ScaleTimesVector(x, betav[0] * EtatoXi, TReta, TReta);
					SharedSpace *SharedTReta = new SharedSpace(TReta);
					etax->AddToTempData("betaTReta", SharedTReta);
				}
			}
			return;
		}
		printf("Warning: The differentiated retraction has not been implemented!\n");
		xix->CopyTo(result);
	};

	double L2Sphere::Beta(Variable *x, Vector *etax) const
	{
		return 1;
	};

	void L2Sphere::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		if (!etax->TempDataExist("xdydn2"))
		{
			Vector *xdy = x->ConstructEmpty();
			SharedSpace *Sharedxdy = new SharedSpace(xdy);
			VectorAddVector(x, x, y, xdy);
			ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
			etax->AddToTempData("xdydn2", Sharedxdy);
		}
		const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
		Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
		scalarVectorAddVector(x, -2.0 * Metric(x, xix, y), xdydn2, xix, result);
	};

	void L2Sphere::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (!etax->TempDataExist("xdydn2"))
		{
			Vector *xdy = x->ConstructEmpty();
			SharedSpace *Sharedxdy = new SharedSpace(xdy);
			VectorAddVector(x, x, y, xdy);
			ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
			etax->AddToTempData("xdydn2", Sharedxdy);
		}
		const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
		Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
		scalarVectorAddVector(x, -2.0 * Metric(x, xiy, x), xdydn2, xiy, result);
	};

	void L2Sphere::HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (!etax->TempDataExist("xdydn2"))
		{
			Vector *xdy = x->ConstructEmpty();
			SharedSpace *Sharedxdy = new SharedSpace(xdy);
			VectorAddVector(x, x, y, xdy);
			ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
			etax->AddToTempData("xdydn2", Sharedxdy);
		}
		const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
		Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
		const double *xdydn2TV = xdydn2->ObtainReadData();

		integer ell = Hx->Getsize()[0];
		integer length = etax->Getlength();
		const double *M = Hx->ObtainReadData();
		double *Hxpy = new double[ell];

		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		integer inc = 1, N = ell;
		// Hxpy <- M(:, start : start + length - 1) * xdydn2TV, 
		// details: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
		dgemv_(transn, &N, &length, &one, const_cast<double *> (M + start * N), &N, const_cast<double *> (xdydn2TV), &inc, &zero, Hxpy, &inc);

		double scalar = -2.0;
		Hx->CopyTo(result);

		Variable *xflat = x->ConstructEmpty();
		x->CopyTo(xflat);
		double *xflatptr = xflat->ObtainWritePartialData();
		xflatptr[0] /= (2 * (n - 1));
		xflatptr[n - 1] /= (2 * (n - 1));
		for (integer i = 1; i < n - 1; i++)
		{
			xflatptr[i] /= (n - 1);
		}
		double *resultL = result->ObtainWritePartialData();
		// resultL(:, start : start + length - 1) <- scalar * Hxpy * xflatptr^T + resultL(:, start : start + length - 1)
		// details: http://www.netlib.org/lapack/explore-html/dc/da8/dger_8f.html
		dger_(&length, &N, &scalar, Hxpy, &inc, xflatptr, &inc, resultL + start * N, &N);
		delete[] Hxpy;
		delete xflat;
	};

	void L2Sphere::TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (!etax->TempDataExist("xdydn2"))
		{
			Vector *xdy = x->ConstructEmpty();
			SharedSpace *Sharedxdy = new SharedSpace(xdy);
			VectorAddVector(x, x, y, xdy);
			ScaleTimesVector(x, 1.0 / Metric(x, xdy, xdy), xdy, xdy);
			etax->AddToTempData("xdydn2", Sharedxdy);
		}

		integer ell = Hx->Getsize()[0];
		integer length = etax->Getlength();
		const double *M = Hx->ObtainReadData();
		double *Hty = new double[ell];

		Variable *yflat = y->ConstructEmpty();
		y->CopyTo(yflat);
		double *yflatptr = yflat->ObtainWritePartialData();
		yflatptr[0] /= (2 * (n - 1));
		yflatptr[n - 1] /= (2 * (n - 1));
		for (integer i = 1; i < n - 1; i++)
		{
			yflatptr[i] /= (n - 1);
		}

		char *transt = const_cast<char *> ("t");
		double one = 1, zero = 0;
		integer inc = 1, N = ell;
		// Hty <- M(start : start + length - 1, :)^T * yflatptr
		// details: http://www.netlib.org/lapack/explore-html/dc/da8/dgemv_8f.html
		dgemv_(transt, &length, &N, &one, const_cast<double *> (M + start), &N, yflatptr, &inc, &zero, Hty, &inc);

		double scalar = -2.0;
		Hx->CopyTo(result);


		const SharedSpace *Sharedxdydn2 = etax->ObtainReadTempData("xdydn2");
		Vector *xdydn2 = Sharedxdydn2->GetSharedElement();
		const double *xdydn2TV = xdydn2->ObtainReadData();

		double *resultL = result->ObtainWritePartialData();

		// resultL(start : start + length - 1, :) <- scalar * xdydn2TV * Hty^T + resultL(start : start + length - 1, :)
		// details: http://www.netlib.org/lapack/explore-html/dc/da8/dger_8f.html
		dger_(&length, &N, &scalar, const_cast<double *> (xdydn2TV), &inc, Hty, &inc, resultL + start, &N);
		delete[] Hty;
		delete yflat;
	};

	void L2Sphere::TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
	{
		HInvTran(x, etax, y, Hx, 0, etax->Getlength(), result);
		TranH(x, etax, y, result, 0, etax->Getlength(), result);
	};

	void L2Sphere::ObtainEtaxFlat(Variable *x, Vector *etax, Vector *etaxflat) const
	{
		etax->CopyTo(etaxflat);
		double *etaxflatTV = etaxflat->ObtainWritePartialData();
		double intv = 1.0 / (n - 1);
		ScaleTimesVector(x, intv, etaxflat, etaxflat);
		etaxflatTV[0] /= 2;
		etaxflatTV[n - 1] /= 2;
	};

	void L2Sphere::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		printf("Routine of obtaining intrinsic representations has not been done!\n");
	};

	void L2Sphere::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		printf("Routine of obtaining extrinsic representations has not been done!\n");
	};

	void L2Sphere::IntrProjection(Variable *x, Vector *v, Vector *result) const
	{
		v->CopyTo(result);
	};

	void L2Sphere::ExtrProjection(Variable *x, Vector *v, Vector *result) const
	{
		const double *xl = x->ObtainReadData();
		double nume = Metric(x, x, v);
		scalarVectorAddVector(x, -nume, x, v, result);
	};

	void L2Sphere::CheckParams(void) const
	{
		std::string Repa2NSMetricnames[L2SPHEREMETRICLENGTH] = { "TRAPEZOID" };
		std::string Repa2NSRetractionnames[L2SPHERERETRACTIONLENGTH] = { "NORMALIZED" };
		std::string Repa2NSVectorTransportnames[L2SPHEREVECTORTRANSPORTLENGTH] = { "L2SPARALLELTRANSLATION" };
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("n             :%15d,\t", n);
		printf("metric        :%15s\n", Repa2NSMetricnames[metric].c_str());
		printf("retraction    :%15s,\t", Repa2NSRetractionnames[retraction].c_str());
		printf("VecTran       :%15s\n", Repa2NSVectorTransportnames[VecTran].c_str());
	};

	void L2Sphere::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be required.
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		ExtrProjection(x, egf, gf);
	};

	void L2Sphere::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		const double *xptr = x->ObtainReadData();
		Variable *xcubed = x->ConstructEmpty();
		SharedSpace *Sharedxcubed = new SharedSpace(xcubed);
		double *xcubedptr = xcubed->ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
		{
			xcubedptr[i] = xptr[i] * xptr[i] * xptr[i];
		}
		double a1 = Metric(x, xcubed, xcubed);

		const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
		Vector *egfVec = Sharedegf->GetSharedElement();
		double a2 = Metric(x, egfVec, xcubed);

		Vector *x2etax = etax->ConstructEmpty();
		double *x2etaxptr = x2etax->ObtainWriteEntireData();
		const double *etaxptr = etax->ObtainReadData();
		for (integer i = 0; i < n; i++)
		{
			x2etaxptr[i] = xptr[i] * xptr[i] * etaxptr[i];
		}
		scalarVectorAddVector(x, -3.0 * a2 / a1, x2etax, exix, xix);
		delete x2etax;
		ExtrProjection(x, xix, xix);
	};
}; /*end of ROPTLIB namespace*/
