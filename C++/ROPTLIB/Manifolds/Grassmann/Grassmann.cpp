
#include "Manifolds/Grassmann/Grassmann.h"

/*Define the namespace*/
namespace ROPTLIB{

	//==============================================================================================
	Grassmann::Grassmann(integer inn, integer inp)
	{
		// public parameter
		HasHHR = false;

		// Only some combinations exist.
		IsIntrApproach = true;
		UpdBetaAlone = false;

		// Status of locking condition
		HasLockCon = false;

		// Fixed parameters
		n = inn;
		p = inp;
		ExtrinsicDim = n * p;
		IntrinsicDim = (n - p) * p;//--- n * p - p * (p + 1) / 2;
		name.assign("Grassmann");
		EMPTYEXTR = new GrassVector(n, p);
		EMPTYINTR = new GrassVector(IntrinsicDim);
	};

	Grassmann::~Grassmann(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	void Grassmann::CheckParams(void) const
	{
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("n             :%15d,\t", n);
		printf("p             :%15d\n", p);
	};

	void Grassmann::ExtrProjection(Variable *x, Vector *v, Vector *result) const
	{
		integer N = n, P = p, inc = 1, Length = N * P;
		double *UtV = new double[P * P];
		const double *U = x->ObtainReadData();
		const double *V = v->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		// UtV = U^T * V, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &P, &P, &N, &GLOBAL::DONE, const_cast<double *> (U), &N, const_cast<double *> (V), &N, &GLOBAL::DZERO, UtV, &P);

		// resultTV <- V, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		if (V != resultTV)
			dcopy_(&Length, const_cast<double *> (V), &inc, resultTV, &GLOBAL::IONE);
		double negone = -1;
		// resultTV = resultTV - U * UtV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &P, &GLOBAL::DNONE, const_cast<double *> (U), &N, UtV, &P, &GLOBAL::DONE, resultTV, &N);
		delete[] UtV;
	};

	void Grassmann::Retraction(Variable *x, Vector *etax, Variable *result) const
	{ // qf retraction: result = R_x(etax) = qf(x + etax)
		const double *U = x->ObtainReadData();
		const double *V;
		Vector *exetax = nullptr;
		if (IsIntrApproach)
		{
			exetax = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(x, etax, exetax);
			V = exetax->ObtainReadData();
		}
		else
		{
			V = etax->ObtainReadData();
		}
		double *resultM = result->ObtainWriteEntireData();
		SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
		double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
		SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1]);
		double *tau = HHRTau->ObtainWriteEntireData();

		integer N = x->Getsize()[0], P = x->Getsize()[1], Length = N * P, inc = 1;
		double one = 1, zero = 0;
		// ptrHHR <- V, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (V), &inc, ptrHHR, &inc);
		// ptrHHR <- U + ptrHHR, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&Length, &one, const_cast<double *> (U), &inc, ptrHHR, &inc);
		integer *jpvt = new integer[P];
		integer info;
		integer lwork = -1;
		double lworkopt;
		for (integer i = 0; i < P; i++)
			jpvt[i] = i + 1;
		// compute the space required in the dgeqp3
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
		// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
		dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
		if (info < 0)
			printf("Error in qr decomposition!\n");
		for (integer i = 0; i < P; i++)
		{
			if (jpvt[i] != (i + 1))
				printf("Error in qf retraction!\n");
		}
		double *signs = new double[P];
		for (integer i = 0; i < P; i++)
			signs[i] = (ptrHHR[i + i * N] >= 0) ? 1 : -1;
		// resultM <- ptrHHR, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, ptrHHR, &inc, resultM, &inc);
		// Generate an orthonormal matrix by using the Householder refections in resultM, output is stored in resultM,
		// details: http://www.netlib.org/lapack/explore-html/d9/d1d/dorgqr_8f.html
		dorgqr_(&N, &P, &P, resultM, &N, tau, work, &lwork, &info);
		if (info < 0)
			printf("Error in forming Q matrix!\n");
		for (integer i = 0; i < P; i++)
			// resultM(:, i) <- signs(i) * resultM(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, signs + i, resultM + i * N, &inc);
		result->AddToTempData("HHR", HouseHolderResult); /*HHR and HHRTau are used to representation (x x_\perp); This orthonormal matrix is
														 used to compute intrinsic representation of tangent vector at x.*/
		result->AddToTempData("HHRTau", HHRTau);
		delete[] jpvt;
		delete[] work;
		delete[] signs;
		if (exetax != nullptr)
			delete exetax;
	};

	void Grassmann::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		Vector *extempx = EMPTYEXTR->ConstructEmpty();
		const double *extempxTV;
		if (IsIntrApproach)
		{
			ObtainExtr(x, xix, extempx);
			extempxTV = extempx->ObtainReadData();
		}
		else
		{
			xix->CopyTo(extempx);
			extempxTV = extempx->ObtainWritePartialData();
		}
		const double *yM = y->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();
		const SharedSpace *HHR = y->ObtainReadTempData("HHR");
		const double *ptrHHR = HHR->ObtainReadData();
		double *YtVRinv = new double[p * p];
		integer inc = 1, N = n, P = p;
		char *left = const_cast<char *> ("r"), *up = const_cast<char *> ("u"),
			*transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t"), *nonunit = const_cast<char *> ("n");
		double one = 1, zero = 0;
		// solving linear system X ptrHHR = extempxTV for X, X is stored in extempxTV,
		// details: http://www.netlib.org/lapack/explore-html/de/da7/dtrsm_8f.html
		dtrsm_(left, up, transn, nonunit, &N, &P, &one, const_cast<double *> (ptrHHR), &N, const_cast<double *> (extempxTV), &N);

		double sign;
		for (integer i = 0; i < P; i++)
		{
			sign = (ptrHHR[i + i * N] >= 0) ? 1 : -1;
			// extempxTV(:, i) <- sign * extempxTV(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &sign, const_cast<double *> (extempxTV + i * N), &inc);
		}
		// YtVRinv <- yM^T * extempxTV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (yM), &N, const_cast<double *> (extempxTV), &N, &zero, YtVRinv, &P);
		for (integer i = 0; i < p; i++)
		{
			YtVRinv[i + p * i] = -YtVRinv[i + p * i];
			for (integer j = i + 1; j < p; j++)
			{
				YtVRinv[i + p * j] = -YtVRinv[j + p * i] - YtVRinv[i + p * j];
				YtVRinv[j + p * i] = 0;
			}
		}
		// extempxTV <- extempxTV + yM * YtVRinv, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &P, &one, const_cast<double *> (yM), &N, YtVRinv, &P, &one, const_cast<double *> (extempxTV), &N);
		if (IsIntrApproach)
		{
			ObtainIntr(y, extempx, result);
		}
		else
		{
			extempx->CopyTo(result);
		}
		delete[] YtVRinv;
		delete extempx;

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
	};

	void Grassmann::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		const double *yM = y->ObtainReadData();
		Vector *exresult = EMPTYEXTR->ConstructEmpty();
		double *exresultTV = exresult->ObtainWriteEntireData();

		Vector *extempy = nullptr;
		const double *extempyTV;
		if (IsIntrApproach)
		{
			extempy = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(y, xiy, extempy);
			extempyTV = extempy->ObtainReadData();
		}
		else
		{
			extempyTV = xiy->ObtainReadData();
		}
		double *ytxiy = new double[p * p];

		char *transt = const_cast<char *> ("t"), *transn = const_cast<char *> ("n");
		integer N = n, P = p, inc = 1;
		double one = 1, zero = 0;
		// ytxiy = yM^T * extempyTV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (yM), &N, const_cast<double *> (extempyTV), &N, &zero, ytxiy, &P);

		// exresultTV = yM * ytxiy, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &P, &one, const_cast<double *> (yM), &N, ytxiy, &P, &zero, exresultTV, &N);
		integer Length = N * P;
		// exresultTV <- extempyTV + exresultTV, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&Length, &one, const_cast<double *> (extempyTV), &inc, exresultTV, &inc);

		const SharedSpace *HHR = y->ObtainReadTempData("HHR");
		const double *ptrHHR = HHR->ObtainReadData();
		double sign;
		for (integer i = 0; i < P; i++)
		{
			sign = (ptrHHR[i + i * N] >= 0) ? 1 : -1;
			// exresultTV(:, i) <- sign * exresultTV(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &sign, exresultTV + i * N, &inc);
		}

		char *left = const_cast<char *> ("r"), *up = const_cast<char *> ("u"), *nonunit = const_cast<char *> ("n");
		// solving linear system X ptrHHR^T = exresultTV for X, X is stored in exresultTV,
		// details: http://www.netlib.org/lapack/explore-html/de/da7/dtrsm_8f.html
		dtrsm_(left, up, transt, nonunit, &N, &P, &one, const_cast<double *> (ptrHHR), &N, exresultTV, &N);

		ExtrProjection(x, exresult, exresult);
		if (IsIntrApproach)
		{
			ObtainIntr(x, exresult, result);
		}
		else
		{
			exresult->CopyTo(result);
		}

		delete[] ytxiy;
		delete exresult;
		if (extempy != nullptr)
			delete extempy;
	};

	double Grassmann::Beta(Variable *x, Vector *etax) const
	{
		if (!HasHHR && !UpdBetaAlone)
			return 1;

		if (!etax->TempDataExist("beta"))
		{ /*In case that beta is not computed, then compute it.*/
			Variable *y = x->ConstructEmpty();
			Vector *xiy = etax->ConstructEmpty();
			Retraction(x, etax, y);
			DiffRetraction(x, etax, y, etax, xiy, true);
			delete y;
			delete xiy;
		}

		/*If the beta has been computed in differentiated retraction, then obtain it.
		Beta should be almost always computed before.*/
		const SharedSpace *beta = etax->ObtainReadTempData("beta");
		const double *betav = beta->ObtainReadData();
		return betav[0];
	};

	void Grassmann::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
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

	void Grassmann::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector *xix, const Problem *prob) const
	{
		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		double one = 1, zero = 0;
		integer inc = 1, N = n, P = p, Length = N * P;
		double *xtegfptr;
		SharedSpace *xtegf;
		if (x->TempDataExist("xtegf"))
		{
			xtegf = const_cast<SharedSpace *> (x->ObtainReadTempData("xtegf"));
			xtegfptr = const_cast<double *> (xtegf->ObtainReadData());
		}
		else
		{
			const double *xxM = x->ObtainReadData();
			const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
			Vector *egfVec = Sharedegf->GetSharedElement();
			const double *egf = egfVec->ObtainReadData();
			xtegf = new SharedSpace(2, p, p);
			xtegfptr = xtegf->ObtainWriteEntireData();
			// xtegfptr <- xxM^T * egf, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
			dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (xxM), &N, const_cast<double *> (egf), &N, &zero, xtegfptr, &P);
		}

		exix->CopyTo(xix);
		double *resultTV = xix->ObtainWritePartialData();
		const double *etaxTV = etax->ObtainReadData();

		double negone = -1;
		// resultTV <- resultTV - etaxTV * xtegfptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (etaxTV), &N, xtegfptr, &P, &one, resultTV, &N);
		ExtrProjection(x, xix, xix);
		if (!x->TempDataExist("xtegf"))
		{
			x->AddToTempData("xtegf", xtegf);
		}
	};

	void Grassmann::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			const double *xM = x->ObtainReadData();
			SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
			double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
			SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1]);
			double *tau = HHRTau->ObtainWriteEntireData();

			integer N = x->Getsize()[0], P = x->Getsize()[1], Length = N * P, inc = 1;
			double one = 1, zero = 0;
			// ptrHHR <- xM, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
			dcopy_(&Length, const_cast<double *> (xM), &inc, ptrHHR, &inc);
			integer *jpvt = new integer[P];
			integer info;
			integer lwork = -1;
			double lworkopt;
			// compute the size of space required in the dgeqp3
			dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
			lwork = static_cast<integer> (lworkopt);
			double *work = new double[lwork];
			for (integer i = 0; i < P; i++)
				jpvt[i] = i + 1;
			// QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
			// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
			dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
			x->AddToTempData("HHR", HouseHolderResult);
			x->AddToTempData("HHRTau", HHRTau);
			if (info < 0)
				printf("Error in qr decomposition!\n");
			for (integer i = 0; i < P; i++)
			{
				if (jpvt[i] != (i + 1))
					printf("Error in qf retraction!\n");
			}
			delete[] jpvt;
			delete[] work;
		}

		const double *xM = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		double *resultTV = result->ObtainWriteEntireData();
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();

		char *transt = const_cast<char *> ("t"), *sidel = const_cast<char *> ("l");
		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		integer info;
		integer lwork = -1;
		double lworkopt;
		double *tempspace = new double[n * p];
		// compute the size of space required in the dormqr
		dormqr_(sidel, transt, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), tempspace, &N, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// tempspace <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&Length, const_cast<double *> (etaxTV), &inc, tempspace, &inc);
		// tempspace <- Q^T * tempspace, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
		dormqr_(sidel, transt, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), tempspace, &N, work, &lwork, &info);

		for (integer i = 0; i < p; i++)
		{
			integer nmp = n - p;
			dcopy_(&nmp, tempspace + p + n* i, &GLOBAL::IONE, resultTV + nmp * i, &GLOBAL::IONE);
		}

		delete[] work;
		delete[] tempspace;
	};

	void Grassmann::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("HHR"))
		{
			const double *xM = x->ObtainReadData();
			SharedSpace *HouseHolderResult = new SharedSpace(2, x->Getsize()[0], x->Getsize()[1]);
			double *ptrHHR = HouseHolderResult->ObtainWriteEntireData();
			SharedSpace *HHRTau = new SharedSpace(1, x->Getsize()[1]);
			double *tau = HHRTau->ObtainWriteEntireData();

			integer N = x->Getsize()[0], P = x->Getsize()[1], Length = N * P, inc = 1;
			double one = 1, zero = 0;
			// ptrHHR <- xM, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
			dcopy_(&Length, const_cast<double *> (xM), &inc, ptrHHR, &inc);
			integer *jpvt = new integer[P];
			integer info;
			integer lwork = -1;
			double lworkopt;
			for (integer i = 0; i < P; i++)
				jpvt[i] = i + 1;
			// compute the size of space required in the dgeqp3
			dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, &lworkopt, &lwork, &info);
			lwork = static_cast<integer> (lworkopt);
			double *work = new double[lwork];
			// QR decomposition for ptrHHR using Householder reflections. Householder reflectors and R are stored in ptrHHR.
			// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
			dgeqp3_(&N, &P, ptrHHR, &N, jpvt, tau, work, &lwork, &info);
			x->AddToTempData("HHR", HouseHolderResult);
			x->AddToTempData("HHRTau", HHRTau);
			if (info < 0)
				printf("Error in qr decomposition!\n");
			for (integer i = 0; i < P; i++)
			{
				if (jpvt[i] != (i + 1))
					printf("Error in qf retraction!\n");
			}
			delete[] jpvt;
			delete[] work;
		}

		const double *xM = x->ObtainReadData();
		const SharedSpace *HHR = x->ObtainReadTempData("HHR");
		const SharedSpace *HHRTau = x->ObtainReadTempData("HHRTau");
		const double *ptrHHR = HHR->ObtainReadData();
		const double *ptrHHRTau = HHRTau->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n"), *sidel = const_cast<char *> ("l");
		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		integer info;
		integer idx = 0;

		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				resultTV[j + i * n] = 0;
			}
			integer nmp = n - p;
			dcopy_(&nmp, const_cast<double *> (intretaxTV)+nmp * i, &GLOBAL::IONE, resultTV + p + n* i, &GLOBAL::IONE);
		}

		double sign;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			// result(i, :) <- sign * result(i, :), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&P, &sign, resultTV + i, &N);
		}
		integer lwork = -1;
		double lworkopt;
		// compute the size of space required in the dormqr
		dormqr_(sidel, transn, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), resultTV, &N, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);

		double *work = new double[lwork];
		// resultTV <- Q * resultTV, where Q is the orthogonal matrix defined as the product elementary reflectors defined by ptrHHR and ptrHHRTau,
		// details: http://www.netlib.org/lapack/explore-html/da/d82/dormqr_8f.html
		dormqr_(sidel, transn, &N, &P, &P, const_cast<double *> (ptrHHR), &N, const_cast<double *> (ptrHHRTau), resultTV, &N, work, &lwork, &info);
		delete[] work;
	};
}; /*end of ROPTLIB namespace*/
