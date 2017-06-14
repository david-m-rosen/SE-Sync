
#include "Manifolds/Stiefel/Stiefel.h"

/*Define the namespace*/
namespace ROPTLIB{

	//==============================================================================================
	Stiefel::Stiefel(integer inn, integer inp)
	{
		// public parameter
		HasHHR = false;

		// Only some combinations exist.
		metric = EUCLIDEAN;
		retraction = QF;
		VecTran = PARALLELIZATION;
		IsIntrApproach = true;
		UpdBetaAlone = false;

		// Status of locking condition
		HasLockCon = false;

		// Fixed parameters
		n = inn;
		p = inp;
		ExtrinsicDim = n * p;
		IntrinsicDim = n * p - p * (p + 1) / 2;
		name.assign("Stiefel");
		EMPTYEXTR = new StieVector(n, p);
		EMPTYINTR = new StieVector(IntrinsicDim);
	};

	Stiefel::~Stiefel(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	// Choose the default parameters
	void Stiefel::ChooseStieParamsSet1(void)
	{
		metric = EUCLIDEAN;
		retraction = QF;
		VecTran = PARALLELIZATION;
		IsIntrApproach = true;
		UpdBetaAlone = false;
		HasHHR = false;
		HasLockCon = false;
	};

	// Choose the parameters
	void Stiefel::ChooseStieParamsSet2(void)
	{
		metric = EUCLIDEAN;
		retraction = CONSTRUCTED;
		VecTran = PARALLELIZATION;
		IsIntrApproach = true;
		UpdBetaAlone = false;
		HasHHR = false;
		HasLockCon = true;
	};

	void Stiefel::ChooseStieParamsSet3(void)
	{
		metric = EUCLIDEAN;
		retraction = QF;
		VecTran = PROJECTION;
		IsIntrApproach = false;
		UpdBetaAlone = false;
		HasHHR = false;
		HasLockCon = false;
	};

	void Stiefel::ChooseStieParamsSet4(void)
	{
		metric = EUCLIDEAN;
		retraction = CAYLEYR;
		VecTran = CAYLEYVT;
		IsIntrApproach = false;
		UpdBetaAlone = false;
		HasHHR = false;
		HasLockCon = false;
	};

	void Stiefel::CheckParams(void) const
	{
		std::string StieMetricnames[STIEMETRICLENGTH] = { "EUCLIDEAN", "CANONICAL" };
		std::string StieRetractionnames[STIERETRACTIONLENGTH] = { "QF", "POLAR", "EXP", "CONSTRUCTED", "CAYLEYR", "PROXSTIE" };
		std::string StieVectorTransportnames[STIEVECTORTRANSPORTLENGTH] = { "PARALLELIZATION", "RIGGING", "PARALLELTRANSLATION", "PROJECTION", "CAYLEYVT" };
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("n             :%15d,\t", n);
		printf("p             :%15d\n", p);
		printf("metric        :%15s,\t", StieMetricnames[metric].c_str());
		printf("retraction    :%15s\n", StieRetractionnames[retraction].c_str());
		printf("VecTran       :%15s\n", StieVectorTransportnames[VecTran].c_str());
	};

	void Stiefel::IntrProjection(Variable *x, Vector *etax, Vector *result) const
	{
		etax->CopyTo(result);
	};

	void Stiefel::ExtrProjection(Variable *x, Vector *v, Vector *result) const
	{
		integer N = n, P = p, inc = 1, Length = N * P;
		double *symUtV = new double[P * P];
		const double *U = x->ObtainReadData();
		const double *V = v->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		double one = 1, zero = 0;
		// symUtV = U^T * V, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (U), &N, const_cast<double *> (V), &N, &zero, symUtV, &P);
		for (integer i = 0; i < P; i++)
		{
			for (integer j = i + 1; j < P; j++)
			{
				symUtV[i + j * P] += symUtV[j + i * P];
				symUtV[i + j * P] /= 2.0;
				symUtV[j + i * P] = symUtV[i + j * P];
			}
		}
		// resultTV <- V, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		if (V != resultTV)
			dcopy_(&Length, const_cast<double *> (V), &inc, resultTV, &inc);
		double negone = -1;
		// resultTV = resultTV - U * symUtV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (U), &N, symUtV, &P, &one, resultTV, &N);
		delete[] symUtV;
	};

	//=================================================================================
	double Stiefel::Metric(Variable *x, Vector *etax, Vector *xix) const
	{
		if (metric == EUCLIDEAN)
			return Manifold::Metric(x, etax, xix);
		printf("Error: Metric has not been done!\n");
		return 0;
	};

	void Stiefel::Projection(Variable *x, Vector *v, Vector *result) const
	{
		if (IsIntrApproach)
			IntrProjection(x, v, result);
		else
			ExtrProjection(x, v, result);
	};

	void Stiefel::Retraction(Variable *x, Vector *etax, Variable *result) const
	{
		if (retraction == QF)
			return qfRetraction(x, etax, result);

		if (retraction == CONSTRUCTED)
			return ConRetraction(x, etax, result);

		if (retraction == CAYLEYR)
			return CayleyRetraction(x, etax, result);

		printf("Error: Retraction has not been done!\n");
	};

	void Stiefel::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (retraction == QF)
			return qfcoTangentVector(x, etax, y, xiy, result);

		if (retraction == CONSTRUCTED)
			return ConcoTangentVector(x, etax, y, xiy, result);

		if (retraction == CAYLEYR)
			return CayleycoTangentVector(x, etax, y, xiy, result);

		printf("Error: coTangentVector has not been done!\n");
	};

	void Stiefel::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (retraction == QF)
			return DiffqfRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

		if (retraction == CONSTRUCTED)
			return DiffConRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

		if (retraction == CAYLEYR)
			return DiffCayleyRetraction(x, etax, y, xix, result, IsEtaXiSameDir);

		printf("Error: DiffRetraction has not been done!\n");
	};

	double Stiefel::Beta(Variable *x, Vector *etax) const
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

	void Stiefel::VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		if (VecTran == PARALLELIZATION && !HasHHR)
		{
			return Manifold::VectorTransport(x, etax, y, xix, result);
		}

		if (VecTran == PROJECTION && !HasHHR)
			return ExtrProjection(y, xix, result);

		if (VecTran == CAYLEYVT && !HasHHR)
			return CayleyVectorTransport(x, etax, y, xix, result);

		if (HasHHR)
			return LCVectorTransport(x, etax, y, xix, result);

		printf("Error: VectorTransport has not been done!\n");
	};

	void Stiefel::InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		if (VecTran == PARALLELIZATION && !HasHHR)
			return Manifold::InverseVectorTransport(x, etax, y, xiy, result);

		if (VecTran == PROJECTION && !HasHHR)
		{
			printf("Stiefel::InverseVectorTransport: inverse vector transport by projection has not been done!\n");
			return Manifold::InverseVectorTransport(x, etax, y, xiy, result);
		}

		if (VecTran == CAYLEYVT && !HasHHR)
			return CayleyInverseVectorTransport(x, etax, y, xiy, result);

		if (HasHHR)
			return LCInverseVectorTransport(x, etax, y, xiy, result);

		printf("Error: InverseVectorTransport has not been done!\n");
	};

	void Stiefel::HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == PARALLELIZATION && !HasHHR)
			return Manifold::HInvTran(x, etax, y, Hx, start, end, result);

		if (VecTran == PROJECTION && !HasHHR)
		{
			printf("Stiefel::HInvTran for vector transport by projection has not been done!\n");
			return Manifold::HInvTran(x, etax, y, Hx, start, end, result);
		}

		if (VecTran == CAYLEYVT && !HasHHR)
		{
			printf("Stiefel::HInvTran for Cayley vector transport has not been done!\n");
			return Manifold::HInvTran(x, etax, y, Hx, start, end, result);
		}

		if (HasHHR)
			return LCHInvTran(x, etax, y, Hx, start, end, result);

		printf("Error: HInvTran has not been done!\n");
	};

	void Stiefel::TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const
	{
		if (VecTran == PARALLELIZATION && !HasHHR)
			return Manifold::TranH(x, etax, y, Hx, start, end, result);

		if (VecTran == PROJECTION && !HasHHR)
		{
			printf("Stiefel::TranH for vector transport by projection has not been done!\n");
			return Manifold::TranH(x, etax, y, Hx, start, end, result);
		}

		if (VecTran == CAYLEYVT && !HasHHR)
		{
			printf("Stiefel::TranH for Cayley vector transport has not been done!\n");
			return Manifold::TranH(x, etax, y, Hx, start, end, result);
		}

		if (HasHHR)
			return LCTranH(x, etax, y, Hx, start, end, result);

		printf("Error: TranH has not been done!\n");
	};

	void Stiefel::TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const
	{
		if (VecTran == PARALLELIZATION && !HasHHR)
			return Manifold::TranHInvTran(x, etax, y, Hx, result);

		if (VecTran == PROJECTION && !HasHHR)
		{
			printf("Stiefel::TranHInvTran for vector transport by projection has not been done!\n");
			return Manifold::TranHInvTran(x, etax, y, Hx, result);
		}

		if (VecTran == CAYLEYVT && !HasHHR)
		{
			printf("Stiefel::TranHInvTran for Cayley vector transport has not been done!\n");
			return Manifold::TranHInvTran(x, etax, y, Hx, result);
		}

		if (HasHHR)
			return LCTranHInvTran(x, etax, y, Hx, result);

		printf("Error: TranHInvTran has not been done!\n");
	};

	void Stiefel::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		if (metric == EUCLIDEAN)
		{
			if (prob->GetUseHess())
			{
				Vector *segf = egf->ConstructEmpty();
				segf->NewMemoryOnWrite(); // I don't remember the reason. It seems to be necessary.
				egf->CopyTo(segf);
				SharedSpace *Sharedegf = new SharedSpace(segf);
				x->AddToTempData("EGrad", Sharedegf);
			}
			ExtrProjection(x, egf, gf);
			return;
		}
		printf("Warning:The function converting Eucidean Gradient to Riemannian Gradient has not been done!\n");
	};

	void Stiefel::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector *xix, const Problem *prob) const
	{
		if (metric == EUCLIDEAN)
		{
			char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
			double one = 1, zero = 0;
			integer inc = 1, N = n, P = p, Length = N * P;
			double *symxtegfptr;
			SharedSpace *symxtegf;
			if (x->TempDataExist("symxtegf"))
			{
				symxtegf = const_cast<SharedSpace *> (x->ObtainReadTempData("symxtegf"));
				symxtegfptr = const_cast<double *> (symxtegf->ObtainReadData());
			}
			else
			{
				const double *xxM = x->ObtainReadData();
				const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
				Vector *egfVec = Sharedegf->GetSharedElement();
				const double *egf = egfVec->ObtainReadData();
				symxtegf = new SharedSpace(2, p, p);
				symxtegfptr = symxtegf->ObtainWriteEntireData();
				// symxtegfptr <- xxM^T * egf, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
				dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (xxM), &N, const_cast<double *> (egf), &N, &zero, symxtegfptr, &P);
				for (integer i = 0; i < p; i++)
				{
					for (integer j = i + 1; j < p; j++)
					{
						symxtegfptr[i + j * p] += symxtegfptr[j + i * p];
						symxtegfptr[i + j * p] /= 2.0;
						symxtegfptr[j + i * p] = symxtegfptr[i + j * p];
					}
				}
			}

			exix->CopyTo(xix);
			double *resultTV = xix->ObtainWritePartialData();
			const double *etaxTV = etax->ObtainReadData();

			double negone = -1;
			// resultTV <- resultTV - etaxTV * symxtegfptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
			dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (etaxTV), &N, symxtegfptr, &P, &one, resultTV, &N);
			ExtrProjection(x, xix, xix);
			if (!x->TempDataExist("symxtegf"))
			{
				x->AddToTempData("symxtegf", symxtegf);
			}
			return;
		}
		printf("Warning:The function converting action of Eucidean Hessian to action of Riemannian Hessian has not been done!\n");
	};

	void Stiefel::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		if (retraction == QF || retraction == PROXSTIE)
			ObtainIntrHHR(x, etax, result);
		else
			if (retraction == CONSTRUCTED)
				ObtainIntrSquare(x, etax, result);
			else
				printf("Warning: computing intrinsic representation from extrinsic has not been implemented!\n");
	};

	void Stiefel::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		if (retraction == QF || retraction == PROXSTIE)
			ObtainExtrHHR(x, intretax, result);
		else
			if (retraction == CONSTRUCTED)
				ObtainExtrSquare(x, intretax, result);
			else
				printf("Warning: computing extrinsic representation from intrinsic has not been implemented!\n");
	};

	void Stiefel::qfRetraction(Variable *x, Vector *etax, Variable *result) const
	{
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

		///*for debug*/
		//double *FullOrth = new double[N * N];
		//dcopy_(&Length, ptrHHR, &inc, FullOrth, &inc);
		//dorgqr_(&N, &N, &P, FullOrth, &N, tau, work, &lwork, &info);
		//ForDebug::Print("FullOrth:", FullOrth, N, N);
		//delete[] FullOrth;
		///*end for debug*/

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

	void Stiefel::qfcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
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
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i; j < p; j++)
			{
				ytxiy[i + j * p] = -ytxiy[i + j * p];
			}
		}
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

	void Stiefel::DiffqfRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
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

	void Stiefel::ObtainIntrHHR(Variable *x, Vector *etax, Vector *result) const
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
		double sign;
		for (integer i = 0; i < p; i++)
		{
			sign = (ptrHHR[i + n * i] >= 0) ? 1 : -1;
			dscal_(&P, &sign, tempspace + i, &N);
		}
		double r2 = sqrt(2.0);
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[idx] = r2 * (tempspace[j + i * n] - tempspace[i + j * n]) / 2;
				//resultTV[idx] = r2 * tempspace[j + i * n];
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[idx] = tempspace[j + i * n];
				idx++;
			}
		}

		delete[] work;
		delete[] tempspace;
	};

	void Stiefel::ObtainExtrHHR(Variable *x, Vector *intretax, Vector *result) const
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
		double r2 = sqrt(2.0);
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			resultTV[i + i * n] = 0;
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx] / r2;
				resultTV[i + j * n] = -resultTV[j + i * n];
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx];
				idx++;
			}
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

	void Stiefel::ConRetraction(Variable *x, Vector *etax, Variable *result) const
	{ // only accept intrinsic approach
		const double *V = nullptr;
		Vector *inetax = EMPTYINTR->ConstructEmpty();
		if (IsIntrApproach)
		{
			V = etax->ObtainReadData();
		}
		else
		{
			ObtainIntr(x, etax, inetax);
			V = inetax->ObtainReadData();
		}
		double *M = new double[3 * n * n + 2 * n];
		double *wr = M + n * n;
		double *wi = wr + n;
		double *Vs = wi + n;
		double *VsT = Vs + n * n;

		double r2 = sqrt(2.0);
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			M[i + i * n] = 0;
			for (integer j = i + 1; j < p; j++)
			{
				M[j + i * n] = V[idx] / r2;
				M[i + j * n] = -M[j + i * n];
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				M[j + i * n] = V[idx];
				M[i + j * n] = -V[idx];
				idx++;
			}
		}
		delete inetax;
		for (integer i = p; i < n; i++)
		{
			for (integer j = p; j < n; j++)
			{
				M[j + i * n] = 0;
			}
		}
		char *jobv = const_cast<char *> ("V"), *sortn = const_cast<char *> ("N");
		integer N = n, P = p, NmP = n - p, sdim, info;
		integer lwork = -1;
		double lworkopt;

		// compute the size of space required in the dgees
		dgees_(jobv, sortn, nullptr, &N, M, &N, &sdim, wr, wi, Vs, &N, &lworkopt, &lwork, nullptr, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// schur factorization for M, M = Vs T Vs^T. For output, T is stored in M,
		// details: http://www.netlib.org/lapack/explore-html/d1/d39/dgees_8f.html
		dgees_(jobv, sortn, nullptr, &N, M, &N, &sdim, wr, wi, Vs, &N, work, &lwork, nullptr, &info);

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		double cosv, sinv;
		integer two = 2, inc = 1;
		double one = 1, zero = 0;
		double block[4];
		for (integer i = 0; i < n; i++)
		{
			if (i + 1 < n && fabs(M[i + (i + 1) * n]) > std::numeric_limits<double>::epsilon())
			{
				cosv = cos(M[i + (i + 1) * n]);
				sinv = sin(M[i + (i + 1) * n]);
				block[0] = cosv; block[1] = -sinv; block[2] = sinv; block[3] = cosv;
				// VsT(:, i : i + 1) <- VsT(:, i : i + 1) * block, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
				dgemm_(transn, transn, &N, &two, &two, &one, Vs + i * n, &N, block, &two, &zero, VsT + i * n, &N);
				i++;
			}
			else
			{
				// VsT(:, i) <- Vs(:, i), details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
				dcopy_(&N, Vs + i * n, &inc, VsT + i * n, &inc);
			}
		}
		// M <- VsT * Vs^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &N, &N, &N, &one, VsT, &N, Vs, &N, &zero, M, &N);

		if (!x->TempDataExist("Perp"))
		{
			ObtainPerp(x);
		}
		const SharedSpace *SharedSpacePerp = x->ObtainReadTempData("Perp");
		const double *Perp = SharedSpacePerp->ObtainReadData();

		const double *xM = x->ObtainReadData();
		double *resultM = result->ObtainWriteEntireData();
		SharedSpace *ResultSharedPerp = new SharedSpace(2, n, n - p);
		double *ResultPerp = ResultSharedPerp->ObtainWriteEntireData();

		// resultM(0 : p-1, 0 : p-1) <- xM(0 : p-1, 0 : p-1) * M(0 : p-1, 0 : p-1), 
		// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &P, &P, &one, const_cast<double *> (xM), &N, M, &N, &zero, resultM, &N);
		// resultM(0 : p-1, 0 : p-1) <- resultM(0 : p-1, 0 : p-1) + Perp(0 : p-1, 0:n-p-1) * M(p : N - 1, 0 : p-1), 
		// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &P, &NmP, &one, const_cast<double *> (Perp), &N, M + p, &N, &one, resultM, &N);

		// resultM(p : N - 1, 0 : p-1) <- xM(p : N - 1, 0:p-1) * M(0:p-1, 0:p-1), 
		// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &NmP, &P, &P, &one, const_cast<double *> (xM + p), &N, M, &N, &zero, resultM + p, &N);
		// resultM(p : N - 1, 0 : p-1) = resultM(p : N - 1, 0 : p-1) + Perp(p : N-1, 0:n-p-1) * M(p : N-1, 0 : p-1), 
		// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &NmP, &P, &NmP, &one, const_cast<double *> (Perp + p), &N, M + p, &N, &one, resultM + p, &N);

		// ResultPerp(0:p-1, 0:n-p-1) <- xM(0:p-1, 0:p-1) * M(0:p-1, p : N-1), 
		// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &NmP, &P, &one, const_cast<double *> (xM), &N, M + n * p, &N, &zero, ResultPerp, &N);
		// ResultPerp(0:p-1, 0:n-p-1) <- Perp(0:p-1, :) * M(p : n-1, p : n-1), 
		// details about dgemm: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &NmP, &NmP, &one, const_cast<double *> (Perp), &N, M + n * p + p, &N, &one, ResultPerp, &N);

		// ResultPerp(p : n-1, 0:n-p-1) <- xM(p:n-1, 0:p-1) * M(0 : p-1, p : n-1), 
		// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &NmP, &NmP, &P, &one, const_cast<double *> (xM + p), &N, M + n * p, &N, &zero, ResultPerp + p, &N);
		// ResultPerp(p : n-1, 0:n-p-1) <- Perp(p:n-1, 0:n-p-1) * M(p : n-1, p : n-1), 
		// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &NmP, &NmP, &NmP, &one, const_cast<double *> (Perp + p), &N, M + n * p + p, &N, &one, ResultPerp + p, &N);

		result->AddToTempData("Perp", ResultSharedPerp);

		delete[] work;
		delete[] M;
	};

	void Stiefel::ConcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		xiy->CopyTo(result);
		printf("The cotangent vector for the constructed retraction has not been implemented!\n");
	};

	void Stiefel::DiffConRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsEtaXiSameDir)
		{
			Manifold::VectorTransport(x, etax, y, xix, result);

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
		printf("Warning: The differentiated retraction of the constructed retraction has not been implemented!\n");
		xix->CopyTo(result);
	};

	void Stiefel::CayleyRetraction(Variable *x, Vector *etax, Variable *result) const
	{//assume extrinsic representation is used for etax
		const double *xptr = x->ObtainReadData();
		const double *etaxptr = etax->ObtainReadData();
		SharedSpace *SharedSpaceMk2 = new SharedSpace(2, 2 * p, 2 * p);
		SharedSpace *SharedSpaceMk3 = new SharedSpace(2, 2 * p, p);
		SharedSpace *SharedSpaceLUP = new SharedSpace(1, 4 * p * p + 2 * p);
		SharedSpace *SharedSpaceU = new SharedSpace(2, n, 2 * p);

		double *Mk2ptr = SharedSpaceMk2->ObtainWriteEntireData();
		double *Mk3ptr = SharedSpaceMk3->ObtainWriteEntireData();
		double *LUPptr = SharedSpaceLUP->ObtainWriteEntireData();
		double *Uptr = SharedSpaceU->ObtainWriteEntireData();

		integer P = p, N = n, P2 = p * 2;
		double half = 0.5;
		/*1-1 block of Mk2*/
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &half, const_cast<double *> (xptr), &N, const_cast<double *> (etaxptr), &N, &GLOBAL::DZERO, Mk2ptr, &P2);
		/*1-2 block of Mk2*/
		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < P2; j++)
			{
				Mk2ptr[i + j * P2] = 0;
			}
			Mk2ptr[i + (i + p) * P2] = 1;
		}
		/*2-1 block of Mk2*/
		double three = 3.0;
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &P, &three, Mk2ptr, &P2, Mk2ptr, &P2, &GLOBAL::DZERO, Mk2ptr + p, &P2);
		dgemm_(GLOBAL::T, GLOBAL::N, &P, &P, &N, &GLOBAL::DNONE, const_cast<double *> (etaxptr), &N, const_cast<double *> (etaxptr), &N, &GLOBAL::DONE, Mk2ptr + p, &P2);
		/*2-2 block of Mk2*/
		for (integer i = 0; i < p; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				Mk2ptr[i + p + (j + p) * P2] = -Mk2ptr[j + i * P2];
			}
		}
		double *Mk1ptr = Mk2ptr + P2 * p;
		/*for computing LUP*/
		integer Psquare4 = p * p * 4;
		double nhalf = -0.5;
		for (integer i = 0; i < 4 * p * p; i++)
			LUPptr[i] = 0;
		daxpy_(&Psquare4, &nhalf, Mk2ptr, &GLOBAL::IONE, LUPptr, &GLOBAL::IONE);
		for (integer i = 0; i < P2; i++)
			LUPptr[i + i * P2] += 1;
		integer info;
		integer *Perm = new integer[P2];
		// LU decomposion for LUPpter, LUPpter = P * L * U, L and U are stored in LUPpter, the permutation matrix is in Perm
		// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
		dgetrf_(&P2, &P2, LUPptr, &P2, Perm, &info);
		for (integer i = 0; i < P2; i++)
			LUPptr[4 * p * p + i] = static_cast<double> (Perm[i]);
		/*Compute Mk3*/
		integer length = 2 * p * p;
		dcopy_(&length, Mk1ptr, &GLOBAL::IONE, Mk3ptr, &GLOBAL::IONE);
		/*Solve the linear system*/
		dgetrs_(GLOBAL::N, &P2, &P, LUPptr, &P2, Perm, Mk3ptr, &P2, &info);
		if (info != 0)
			printf("Warning: dgetrs in Stiefel::CayleyRetraction failed!\n");
		delete[] Perm;
		/*Compute U*/
		length = n * p;
		dcopy_(&length, const_cast<double *>(etaxptr), &GLOBAL::IONE, Uptr, &GLOBAL::IONE);
		dcopy_(&length, const_cast<double *>(xptr), &GLOBAL::IONE, Uptr + n * p, &GLOBAL::IONE);
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P, &GLOBAL::DNONE, const_cast<double *> (xptr), &N, Mk2ptr, &P2, &GLOBAL::DONE, Uptr, &N);
		/*compute the result: x + U * Mk3*/
		x->CopyTo(result);
		double *resultptr = result->ObtainWritePartialData();
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, Uptr, &N, Mk3ptr, &P2, &GLOBAL::DONE, resultptr, &N);

		etax->AddToTempData("Mk2", SharedSpaceMk2);
		etax->AddToTempData("Mk3", SharedSpaceMk3);
		etax->AddToTempData("LUP", SharedSpaceLUP);
		etax->AddToTempData("U", SharedSpaceU);
	};

	void Stiefel::CayleycoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		xiy->CopyTo(result);
		printf("The cotangent vector for the Cayley retraction has not been implemented!\n");
	};

	void Stiefel::DiffCayleyRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		if (IsEtaXiSameDir)
		{
			const SharedSpace *SharedSpaceMk2 = etax->ObtainReadTempData("Mk2"); /*2p by 2p*/
			const SharedSpace *SharedSpaceMk3 = etax->ObtainReadTempData("Mk3"); /*2p by p */
			const SharedSpace *SharedSpaceLUP = etax->ObtainReadTempData("LUP"); /*2p by 2p + 2p*/
			const SharedSpace *SharedSpaceU = etax->ObtainReadTempData("U");     /*n by 2p */

			const double *Mk2ptr = SharedSpaceMk2->ObtainReadData();
			const double *Mk3ptr = SharedSpaceMk3->ObtainReadData();
			const double *LUPptr = SharedSpaceLUP->ObtainReadData();
			const double *Uptr = SharedSpaceU->ObtainReadData();

			double *Mk2Mk3 = new double[4 * p * p];
			double *tmp = Mk2Mk3 + 2 * p * p;
			integer P2 = p * 2, P = p, N = n, length = 2 * p * p;
			/*compute Mk2 * Mk3 */
			dgemm_(GLOBAL::N, GLOBAL::N, &P2, &P, &P2, &GLOBAL::DONE, const_cast<double *> (Mk2ptr), &P2, const_cast<double *> (Mk3ptr), &P2, &GLOBAL::DZERO, Mk2Mk3, &P2);
			/*tmp <-- Mk1*/
			dcopy_(&length, const_cast<double *> (Mk2ptr + 2 * p * p), &GLOBAL::IONE, tmp, &GLOBAL::IONE);
			double half = 0.5;
			/*tmp <-- Mk1 + 0.5 * Mk2 * Mk3 */
			daxpy_(&length, &half, Mk2Mk3, &GLOBAL::IONE, tmp, &GLOBAL::IONE);

			/*solve (1 - 0.5Mk2)^{-1} Mk2Mk3*/
			integer info, *perm = new integer[2 * p];
			for (integer i = 0; i < 2 * p; i++)
				perm[i] = static_cast<integer> (LUPptr[4 * p * p + i]);

			/*Solve the linear system*/
			dgetrs_(GLOBAL::N, &P2, &P, const_cast<double *>(LUPptr), &P2, perm, Mk2Mk3, &P2, &info);
			if (info != 0)
				printf("Warning: dgetrs in Stiefel::DiffCayleyRetraction failed!\n");
			delete[] perm;
			/*tmp <--Mk1 + 0.5 * Mk2 * Mk3 + 0.5 (1 - 0.5Mk2)^{-1} Mk2Mk3 */
			daxpy_(&length, &half, Mk2Mk3, &GLOBAL::IONE, tmp, &GLOBAL::IONE);

			double *resultptr = result->ObtainWriteEntireData();
			dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, const_cast<double *> (Uptr), &N, tmp, &P2, &GLOBAL::DZERO, resultptr, &N);

			delete[] Mk2Mk3;
			ScaleTimesVector(x, sqrt(Metric(x, xix, xix) / Metric(x, etax, etax)), result, result);

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
		printf("Warning: The differentiated retraction of the constructed retraction has not been implemented!\n");
		xix->CopyTo(result);
	};

	void Stiefel::CayleyVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const
	{
		const SharedSpace *SharedSpaceLUP = etax->ObtainReadTempData("LUP"); /*2p by 2p + 2p*/
		const SharedSpace *SharedSpaceU = etax->ObtainReadTempData("U");     /*n by 2p */

		const double *LUPptr = SharedSpaceLUP->ObtainReadData();
		const double *Uptr = SharedSpaceU->ObtainReadData();

		double *V = new double[n * 2 * p + 2 * p * p];
		double *Vtxix = V + n * 2 * p;
		integer length = n * p;
		dcopy_(&length, const_cast<double *> (Uptr + n * p), &GLOBAL::IONE, V, &GLOBAL::IONE);
		dcopy_(&length, const_cast<double *> (Uptr), &GLOBAL::IONE, V + n * p, &GLOBAL::IONE);
		dscal_(&length, &GLOBAL::DNONE, V + n * p, &GLOBAL::IONE);
		const double *xixptr = xix->ObtainReadData();
		integer N = n, P = p, P2 = p * 2;
		dgemm_(GLOBAL::T, GLOBAL::N, &P2, &P, &N, &GLOBAL::DONE, V, &N, const_cast<double *> (xixptr), &N, &GLOBAL::DZERO, Vtxix, &P2);

		integer info, *Perm = new integer[P2];
		for (integer i = 0; i < P2; i++)
			Perm[i] = static_cast<integer> (LUPptr[4 * p * p + i]);
		/*Solve the linear system*/
		dgetrs_(GLOBAL::N, &P2, &P, const_cast<double *>(LUPptr), &P2, Perm, Vtxix, &P2, &info);
		if (info != 0)
			printf("Warning: dgetrs in Stiefel::DiffCayleyRetraction failed!\n");
		delete[] Perm;

		xix->CopyTo(result);
		double *resultptr = result->ObtainWritePartialData();
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, const_cast<double *> (Uptr), &N, Vtxix, &P2, &GLOBAL::DONE, resultptr, &N);
		delete[] V;
	};

	void Stiefel::CayleyInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		const SharedSpace *SharedSpaceLUP = etax->ObtainReadTempData("LUP"); /*2p by 2p + 2p*/
		const SharedSpace *SharedSpaceU = etax->ObtainReadTempData("U");     /*n by 2p */

		const double *LUPptr = SharedSpaceLUP->ObtainReadData();
		const double *Uptr = SharedSpaceU->ObtainReadData();

		double *V = new double[n * 2 * p + 2 * p * p];
		double *Utxiy = V + n * 2 * p;
		integer length = n * p;
		dcopy_(&length, const_cast<double *> (Uptr + n * p), &GLOBAL::IONE, V, &GLOBAL::IONE);
		dcopy_(&length, const_cast<double *> (Uptr), &GLOBAL::IONE, V + n * p, &GLOBAL::IONE);
		dscal_(&length, &GLOBAL::DNONE, V + n * p, &GLOBAL::IONE);
		const double *xiyptr = xiy->ObtainReadData();
		integer N = n, P = p, P2 = p * 2;
		dgemm_(GLOBAL::T, GLOBAL::N, &P2, &P, &N, &GLOBAL::DONE, const_cast<double *> (Uptr), &N, const_cast<double *> (xiyptr), &N, &GLOBAL::DZERO, Utxiy, &P2);

		integer info, *Perm = new integer[P2];
		for (integer i = 0; i < P2; i++)
			Perm[i] = static_cast<integer> (LUPptr[4 * p * p + i]);
		/*Solve the linear system*/
		dgetrs_(GLOBAL::T, &P2, &P, const_cast<double *>(LUPptr), &P2, Perm, Utxiy, &P2, &info);
		if (info != 0)
			printf("Warning: dgetrs in Stiefel::DiffCayleyRetraction failed!\n");
		delete[] Perm;

		xiy->CopyTo(result);
		double *resultptr = result->ObtainWritePartialData();
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &P, &P2, &GLOBAL::DONE, V, &N, Utxiy, &P2, &GLOBAL::DONE, resultptr, &N);
		delete[] V;
	};

	void Stiefel::ObtainPerp(Variable *x) const
	{
		const double *xM = x->ObtainReadData();
		SharedSpace *SharedSpacePerp = new SharedSpace(2, n, n - p);
		double *Perp = SharedSpacePerp->ObtainWriteEntireData();
		for (integer i = 0; i < n * (n - p); i++)
			Perp[i] = genrandnormal();
		double *temp = new double[p * (n - p)];
		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		double one = 1, zero = 0, neg_one = -1;
		integer P = p, N = n, NmP = n - p;
		// temp <- xM^T * Perp, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &P, &NmP, &N, &one, const_cast<double *> (xM), &N, Perp, &N, &zero, temp, &P);
		// Perp <- Perp - xM * temp, details about dgemm: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &NmP, &P, &neg_one, const_cast<double *> (xM), &N, temp, &P, &one, Perp, &N);
		delete[] temp;

		integer *jpvt = new integer[NmP];
		integer lwork = 2 * NmP + (1 + NmP) * 64, info; // 64 = INITIALBLOCKSIZE
		double *tau = new double[NmP + lwork];
		double *work = tau + NmP;
		for (integer i = 0; i < NmP; i++)
			jpvt[i] = 0;
		// QR decomposition for matrix Perp
		// details: http://www.netlib.org/lapack/explore-html/db/de5/dgeqp3_8f.html
		dgeqp3_(&N, &NmP, Perp, &N, jpvt, tau, work, &lwork, &info);
		if (info < 0)
			printf("Error in qr decomposition!\n");
		// Generate an orthonormal matrix by using the Householder refections in Perp, output is stored in Perp,
		// details: http://www.netlib.org/lapack/explore-html/d9/d1d/dorgqr_8f.html
		dorgqr_(&N, &NmP, &NmP, Perp, &N, tau, work, &lwork, &info);
		if (info < 0)
			printf("Error in forming Q matrix!\n");
		delete[] jpvt;
		delete[] tau;

		x->AddToTempData("Perp", SharedSpacePerp);
	};

	void Stiefel::ObtainIntrSquare(Variable *x, Vector *etax, Vector *result) const
	{
		if (!x->TempDataExist("Perp"))
		{
			ObtainPerp(x);
		}
		const SharedSpace *SharedSpacePerp = x->ObtainReadTempData("Perp");
		const double *Perp = SharedSpacePerp->ObtainReadData();

		const double *xM = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();
		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		integer N = n, P = p, NmP = n - p;
		double one = 1, zero = 0;
		double *tempspace = new double[n * p];
		// tempspace <- xM^T * etaxTV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (xM), &N, const_cast<double *> (etaxTV), &N, &zero, tempspace, &N);
		// tempspace(p : n-1, :) <- Perp^T * etaxTV, details about dgemm: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &NmP, &P, &N, &one, const_cast<double *> (Perp), &N, const_cast<double *> (etaxTV), &N, &zero, tempspace + p, &N);

		double *resultTV = result->ObtainWriteEntireData();
		double r2 = sqrt(2.0);
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			for (integer j = i + 1; j < p; j++)
			{
				resultTV[idx] = r2 * (tempspace[j + i * n] - tempspace[i + j * n]) / 2;
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				resultTV[idx] = tempspace[j + i * n];
				idx++;
			}
		}

		delete[] tempspace;
	};

	void Stiefel::ObtainExtrSquare(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("Perp"))
		{
			ObtainPerp(x);
		}
		const SharedSpace *SharedSpacePerp = x->ObtainReadTempData("Perp");
		const double *Perp = SharedSpacePerp->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *tempspace = new double[n * p];
		double r2 = sqrt(2.0);
		integer idx = 0;
		for (integer i = 0; i < p; i++)
		{
			tempspace[i + i * n] = 0;
			for (integer j = i + 1; j < p; j++)
			{
				tempspace[j + i * n] = intretaxTV[idx] / r2;
				tempspace[i + j * n] = -tempspace[j + i * n];
				idx++;
			}
		}

		for (integer i = 0; i < p; i++)
		{
			for (integer j = p; j < n; j++)
			{
				tempspace[j + i * n] = intretaxTV[idx];
				idx++;
			}
		}
		double *resultTV = result->ObtainWriteEntireData();
		const double *xM = x->ObtainReadData();

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		integer N = n, P = p, NmP = n - p;
		double one = 1, zero = 0;

		// resultTV <- xM * tempspace, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &P, &one, const_cast<double *> (xM), &N, tempspace, &N, &zero, resultTV, &N);
		// resultTV <- Perp * tempspace(p : n-1, :) + resultTV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &NmP, &one, const_cast<double *> (Perp), &N, tempspace + p, &N, &one, resultTV, &N);

		delete[] tempspace;
	};

	void Stiefel::SetParams(PARAMSMAP params)
	{
		Manifold::SetParams(params);
		PARAMSMAP::iterator iter;
		for (iter = params.begin(); iter != params.end(); iter++)
		{
			if (iter->first == static_cast<std::string> ("ParamSet"))
			{
				switch (static_cast<integer> (iter->second))
				{
				case 1:
					ChooseStieParamsSet1();
					break;
				case 2:
					ChooseStieParamsSet2();
					break;
				default:
					break;
				}
			}
		}
	};
}; /*end of ROPTLIB namespace*/
