
#include "Manifolds/LowRank/LowRank.h"

/*Define the namespace*/
namespace ROPTLIB{

	LowRank::LowRank(integer inm, integer inn, integer inr) : ProductManifold(3,
		new Grassmann(inm, inr), static_cast<integer> (1), new Euclidean(inr, inr), static_cast<integer> (1), new Grassmann(inn, inr), static_cast<integer> (1))
	{
		m = inm;
		n = inn;
		r = inr;
		name.assign("LowRank");
		delete EMPTYEXTR;
		delete EMPTYINTR;
		EMPTYEXTR = new LowRankVector(m, r, r, n, r);
		EMPTYINTR = new LowRankVector(m * r - r * r, 1, r, n * r - r * r, 1);
	};

	LowRank::~LowRank()
	{
		for (integer i = 0; i < numofmani; i++)
		{
			delete manifolds[i];
		}
	};

	double LowRank::ExtrMetric(Variable *x, Vector *etax, Vector *xix) const
	{
		LowRankVector *LRetax = dynamic_cast<LowRankVector *> (etax);
		LowRankVector *LRxix = dynamic_cast<LowRankVector *> (xix);
		LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);

		LowRankVector *cpLRetax = LRetax->ConstructEmpty();
		cpLRetax->NewMemoryOnWrite();
		LRetax->CopyTo(cpLRetax);
		LowRankVector *cpLRxix = LRxix->ConstructEmpty();
		cpLRxix->NewMemoryOnWrite();
		LRxix->CopyTo(cpLRxix);
		const double *dU1 = LRetax->GetElement(0)->ObtainReadData();
		const double *dV1 = LRetax->GetElement(2)->ObtainReadData();
		const double *dU2 = LRxix->GetElement(0)->ObtainReadData();
		const double *dV2 = LRxix->GetElement(2)->ObtainReadData();
		const double *D = LRx->GetElement(1)->ObtainReadData();

		double *cdU1 = cpLRetax->GetElement(0)->ObtainWriteEntireData();
		double *cdV1 = cpLRetax->GetElement(2)->ObtainWriteEntireData();
		double *cdU2 = cpLRxix->GetElement(0)->ObtainWriteEntireData();
		double *cdV2 = cpLRxix->GetElement(2)->ObtainWriteEntireData();
		integer mm = m, nn = n, rr = r;
		dgemm_(GLOBAL::N, GLOBAL::N, &mm, &rr, &rr, &GLOBAL::DONE, const_cast<double *> (dU1), &mm, const_cast<double *> (D), &rr, &GLOBAL::DZERO, cdU1, &mm);
		dgemm_(GLOBAL::N, GLOBAL::N, &mm, &rr, &rr, &GLOBAL::DONE, const_cast<double *> (dU2), &mm, const_cast<double *> (D), &rr, &GLOBAL::DZERO, cdU2, &mm);

		dgemm_(GLOBAL::N, GLOBAL::T, &nn, &rr, &rr, &GLOBAL::DONE, const_cast<double *> (dV1), &nn, const_cast<double *> (D), &rr, &GLOBAL::DZERO, cdV1, &nn);
		dgemm_(GLOBAL::N, GLOBAL::T, &nn, &rr, &rr, &GLOBAL::DONE, const_cast<double *> (dV2), &nn, const_cast<double *> (D), &rr, &GLOBAL::DZERO, cdV2, &nn);
		double result = ProductManifold::Metric(x, cpLRetax, cpLRxix);
		delete cpLRetax;
		delete cpLRxix;
		return result;
	};

	void LowRank::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{
		LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
		LowRankVector *LRetax = dynamic_cast<LowRankVector *> (etax);
		LowRankVector *LRresult = dynamic_cast<LowRankVector *> (result);
		LRresult->NewMemoryOnWrite();

		manifolds[0]->ObtainIntr(LRx->GetElement(0), LRetax->GetElement(0), LRresult->GetElement(0));
		manifolds[1]->ObtainIntr(LRx->GetElement(1), LRetax->GetElement(1), LRresult->GetElement(1));
		manifolds[2]->ObtainIntr(LRx->GetElement(2), LRetax->GetElement(2), LRresult->GetElement(2));

		const double *D = LRx->GetElement(1)->ObtainReadData();
		double *UK = LRresult->GetElement(0)->ObtainWritePartialData();
		double *VK = LRresult->GetElement(2)->ObtainWritePartialData();

		double *UKD = new double[(m - r) * r + (n - r) * r];
		double *VKD = UKD + (m - r) * r;

		integer MmR = m - r, NmR = n - r, R = r, RR = r * r, inc = 1, MmRR = (m - r) * r, NmRR = (n - r) * r;
		double one = 1, zero = 0;
		// UKD <- UK * D, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(GLOBAL::N, GLOBAL::N, &MmR, &R, &R, &one, UK, &MmR, const_cast<double *> (D), &R, &zero, UKD, &MmR);
		// VKD <- VK * D^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(GLOBAL::N, GLOBAL::T, &NmR, &R, &R, &one, VK, &NmR, const_cast<double *> (D), &R, &zero, VKD, &NmR);
		// UK <- UKD, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&MmRR, UKD, &inc, UK, &inc);
		// VK <- VKD, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&NmRR, VKD, &inc, VK, &inc);

		delete[] UKD;
	};

	void LowRank::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
		LowRankVector *LRintretax = dynamic_cast<LowRankVector *> (intretax);
		LowRankVector *LRresult = dynamic_cast<LowRankVector *> (result);
		LRresult->NewMemoryOnWrite();
		LowRankVector *LRintretaxDinv = LRintretax->ConstructEmpty();
		LRintretaxDinv->NewMemoryOnWrite();
		LRintretax->CopyTo(LRintretaxDinv);
		const double *D = LRx->GetElement(1)->ObtainReadData();
		double *UKD = LRintretaxDinv->GetElement(0)->ObtainWritePartialData();
		double *VKD = LRintretaxDinv->GetElement(2)->ObtainWritePartialData();
		double *UK = new double[(m - r) * r + (n - r) * r + r * r];
		double *VK = UK + (m - r) * r;
		double *Dinv = VK + (n - r) * r;

		/*Get LU decomposition of D, which has been done*/
		LUofDinx(x);
		const SharedSpace *SharedSpacetmp = x->ObtainReadTempData("LUofD");
		const double *LUofD = SharedSpacetmp->ObtainReadData();
		integer *IPIV = new integer[r];
		for (integer i = 0; i < r; i++)
			IPIV[i] = static_cast<integer> (LUofD[r*r + i]);

		//- It is not a good approach when r is not much smaller than n or m-------
		integer MmR = m - r, NmR = n - r, R = r, inc = 1, RR = r * r, MmRR = (m - r) * r, NmRR = (n - r) * r, info;
		dcopy_(&RR, const_cast<double *> (LUofD), &inc, Dinv, &inc);
		integer lwork = -1;
		double lworkopt;
		// compute the required space in dgetri
		dgetri_(&R, Dinv, &R, IPIV, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// compute the inverse matrix of Dinv using its LU decomposition information from dgetrf
		// details: http://www.netlib.org/lapack/explore-html/df/da4/dgetri_8f.html
		dgetri_(&R, Dinv, &R, IPIV, work, &lwork, &info);
		delete[] work;
		delete[] IPIV;
		double one = 1, zero = 0;
		// UK = UKD * Dinv, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(GLOBAL::N, GLOBAL::N, &MmR, &R, &R, &one, UKD, &MmR, Dinv, &R, &zero, UK, &MmR);
		// VK = VKD * Dinv^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(GLOBAL::N, GLOBAL::T, &NmR, &R, &R, &one, VKD, &NmR, Dinv, &R, &zero, VK, &NmR);
		// UKD <- UK, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&MmRR, UK, &inc, UKD, &inc);
		// VKD <- VK, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&NmRR, VK, &inc, VKD, &inc);

		manifolds[0]->ObtainExtr(LRx->GetElement(0), LRintretaxDinv->GetElement(0), LRresult->GetElement(0));
		manifolds[1]->ObtainExtr(LRx->GetElement(1), LRintretaxDinv->GetElement(1), LRresult->GetElement(1));
		manifolds[2]->ObtainExtr(LRx->GetElement(2), LRintretaxDinv->GetElement(2), LRresult->GetElement(2));

		delete[] UK;
		delete LRintretaxDinv;
	};

	void LowRank::Retraction(Variable *x, Vector *etax, Variable *result) const
	{
		if (IsIntrApproach)
		{
			Vector *exetax = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(x, etax, exetax);
			SetIsIntrApproach(false);
			ProductManifold::Retraction(x, exetax, result);
			SetIsIntrApproach(true);
			delete exetax;
		}
		else
		{
			ProductManifold::Retraction(x, etax, result);
		}
	};

	void LowRank::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		Vector *exxiy = EMPTYEXTR->ConstructEmpty();
		Vector *exresult = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);
		ObtainExtr(y, xiy, exxiy);

		LowRankVariable *LRy = dynamic_cast<LowRankVariable *> (y);
		LowRankVector *LRexxiy = dynamic_cast<LowRankVector *> (exxiy);

		for (integer i = 0; i < numofmani; i++)
		{
			manifolds[i]->SetIsIntrApproach(false);
		}
		double *exxiyU = LRexxiy->GetElement(0)->ObtainWritePartialData();
		double *exxiyD = LRexxiy->GetElement(1)->ObtainWritePartialData();
		double *exxiyV = LRexxiy->GetElement(2)->ObtainWritePartialData();
		const double *Uy = LRy->GetElement(0)->ObtainReadData();
		const double *Dy = LRy->GetElement(1)->ObtainReadData();
		const double *Vy = LRy->GetElement(2)->ObtainReadData();

		double *UxiDDyt = new double[2 * m * r];
		double *Utemp = UxiDDyt + m * r;

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		integer M = m, R = r, N = n, MR = M * R, NR = N * R, inc = 1;
		double one = 1, zero = 0, negone = -1;
		// Utemp <- exxiyU * Dy, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &M, &R, &R, &one, exxiyU, &M, const_cast<double *> (Dy), &R, &zero, Utemp, &M);
		// exxiyU <- Utemp * Dy^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &M, &R, &R, &one, Utemp, &M, const_cast<double *> (Dy), &R, &zero, exxiyU, &M);
		// Utemp <- Uy * exxiyD, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &M, &R, &R, &one, const_cast<double *> (Uy), &M, exxiyD, &R, &zero, Utemp, &M);
		// UxiDDyt <- Utemp * Dy^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &M, &R, &R, &one, Utemp, &M, const_cast<double *> (Dy), &R, &zero, UxiDDyt, &M);
		// exxiyU <- UxiDDyt + exxiyU, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&MR, &one, UxiDDyt, &inc, exxiyU, &inc);

		delete[] UxiDDyt;

		double *VxiDtDy = new double[2 * n * r];
		double *Vtemp = VxiDtDy + n * r;
		// Compute Vtemp <- exxiyV * Dy^T, details about dgemm: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &N, &R, &R, &one, exxiyV, &N, const_cast<double *> (Dy), &R, &zero, Vtemp, &N);
		// Compute exxiyV <- Vtemp * Dy, details about dgemm: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &R, &R, &one, Vtemp, &N, const_cast<double *> (Dy), &R, &zero, exxiyV, &N);
		// Compute Vtemp <- Vy * exxiyD^T, details about dgemm: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &N, &R, &R, &one, const_cast<double *> (Vy), &N, exxiyD, &R, &zero, Vtemp, &N);
		// Compute VxiDtDy <- Vtemp * Dy, details about dgemm: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &R, &R, &one, Vtemp, &N, const_cast<double *> (Dy), &R, &zero, VxiDtDy, &N);
		// exxiyV <- VxiDtDy + exxiyV, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&NR, &one, VxiDtDy, &inc, exxiyV, &inc);

		delete[] VxiDtDy;

		ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
		ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
		ProdVector *prodexetax = dynamic_cast<ProdVector *> (exetax);
		ProdVector *prodexxiy = dynamic_cast<ProdVector *> (exxiy);
		ProdVector *prodexresult = dynamic_cast<ProdVector *> (exresult);
		manifolds[0]->ExtrProjection(prody->GetElement(0), prodexxiy->GetElement(0), prodexxiy->GetElement(0));
		manifolds[2]->ExtrProjection(prody->GetElement(2), prodexxiy->GetElement(2), prodexxiy->GetElement(2));
		prodexresult->NewMemoryOnWrite();

		ProductManifold::coTangentVector(x, exetax, y, exxiy, exresult);

		ExtrProjectionStiePerp(prodx->GetElement(0), prodexresult->GetElement(0), prodexresult->GetElement(0));
		ExtrProjectionStiePerp(prodx->GetElement(2), prodexresult->GetElement(2), prodexresult->GetElement(2));


		//- It is not a good approach if p is not much smaller than n and m-------
		LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
		const double *Dx = LRx->GetElement(1)->ObtainReadData();
		double *dU = prodexresult->GetElement(0)->ObtainWritePartialData();
		double *dV = prodexresult->GetElement(2)->ObtainWritePartialData();
		integer RR = R * R, info;
		integer *IPIV = new integer[r];
		double *Dinv = new double[r * r + m * r + n * r];
		Utemp = Dinv + r * r;
		Vtemp = Utemp + m * r;
		// Dinv <- Dx, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&RR, const_cast<double *> (Dx), &inc, Dinv, &inc);
		// LU decomposion for Dinv, Dinv = P * L * U, L and U are stored in Dinv, the permutation matrix is in IPIV
		// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
		dgetrf_(&R, &R, Dinv, &R, IPIV, &info);
		integer lwork = -1;
		double lworkopt;
		// compute the required space in dgetri
		dgetri_(&R, Dinv, &R, IPIV, &lworkopt, &lwork, &info);
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		// compute the inverse matrix of Dinv using its LU decomposition information from dgetrf
		// details: http://www.netlib.org/lapack/explore-html/df/da4/dgetri_8f.html
		dgetri_(&R, Dinv, &R, IPIV, work, &lwork, &info);
		delete[] work;
		delete[] IPIV;

		// Utemp <- dU * Dinv^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &M, &R, &R, &one, dU, &M, Dinv, &R, &zero, Utemp, &M);
		// dU <- Utemp * Dinv, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &M, &R, &R, &one, Utemp, &M, Dinv, &R, &zero, dU, &M);

		// Vtemp <- dV * Dinv, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &R, &R, &one, dV, &N, Dinv, &R, &zero, Vtemp, &N);
		// dV <- Vtemp * Dinv^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &N, &R, &R, &one, Vtemp, &N, Dinv, &R, &zero, dV, &N);

		delete[] Dinv;
		ObtainIntr(x, exresult, result);
		for (integer i = 0; i < numofmani; i++)
		{
			manifolds[i]->SetIsIntrApproach(true);
		}
		delete exetax;
		delete exxiy;
		delete exresult;
	};

	void LowRank::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		Vector *exxix = EMPTYEXTR->ConstructEmpty();
		Vector *exresult = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);
		ObtainExtr(x, xix, exxix);

		for (integer i = 0; i < numofmani; i++)
		{
			manifolds[i]->SetIsIntrApproach(false);
		}

		ProdVariable *prodx = dynamic_cast<ProdVariable *> (x);
		ProdVector *prodexetax = dynamic_cast<ProdVector *> (exetax);
		ProdVariable *prody = dynamic_cast<ProdVariable *> (y);
		ProdVector *prodexxix = dynamic_cast<ProdVector *> (exxix);
		ProdVector *prodexresult = dynamic_cast<ProdVector *> (exresult);
		prodexresult->NewMemoryOnWrite();
		manifolds[0]->DiffRetraction(prodx->GetElement(0), prodexetax->GetElement(0), prody->GetElement(0),
			prodexxix->GetElement(0), prodexresult->GetElement(0), IsEtaXiSameDir);
		manifolds[1]->DiffRetraction(prodx->GetElement(1), prodexetax->GetElement(1), prody->GetElement(1),
			prodexxix->GetElement(1), prodexresult->GetElement(1), IsEtaXiSameDir);
		manifolds[2]->DiffRetraction(prodx->GetElement(2), prodexetax->GetElement(2), prody->GetElement(2),
			prodexxix->GetElement(2), prodexresult->GetElement(2), IsEtaXiSameDir);

		ObtainIntr(y, exresult, result);
		for (integer i = 0; i < numofmani; i++)
		{
			manifolds[i]->SetIsIntrApproach(true);
		}
		delete exetax;
		delete exxix;
		delete exresult;

		if (IsEtaXiSameDir)
		{
			const double *etaxTV = etax->ObtainReadData();
			const double *xixTV = xix->ObtainReadData();
			double EtatoXi = sqrt(Metric(x, etax, etax) / Metric(x, xix, xix));
			SharedSpace *beta = new SharedSpace(1, 1);
			double *betav = beta->ObtainWriteEntireData();
			betav[0] = sqrt(Metric(x, etax, etax) / Metric(x, result, result)) / EtatoXi;
			etax->AddToTempData("beta", beta);

			Vector *TReta = result->ConstructEmpty();
			result->CopyTo(TReta);
			ScaleTimesVector(x, betav[0] * EtatoXi, TReta, TReta);
			SharedSpace *SharedTReta = new SharedSpace(TReta);
			etax->AddToTempData("betaTReta", SharedTReta);
		}
	};

	void LowRank::ExtrProjection(Variable *x, Vector *etax, Vector *result) const
	{
		Vector *inetax = EMPTYINTR->ConstructEmpty();
		ObtainIntr(x, etax, inetax);
		ObtainExtr(x, inetax, result);
		delete inetax;
		//etax->CopyTo(result);
	};

	void LowRank::ExtrProjectionStiePerp(Variable *x, Vector *v, Vector *result) const
	{
		integer N = x->Getsize()[0], P = x->Getsize()[1], inc = 1, Length = N * P;
		double *UtV = new double[P * P];
		const double *U = x->ObtainReadData();
		const double *V = v->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		double one = 1, zero = 0;
		// UtV <- U^T * V, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (U), &N, const_cast<double *> (V), &N, &zero, UtV, &P);

		// resultTV <- V, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		if (V != resultTV)
			dcopy_(&Length, const_cast<double *> (V), &inc, resultTV, &inc);
		double negone = -1;
		// resultTV = resultTV - U * UtV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (U), &N, UtV, &P, &one, resultTV, &N);
		delete[] UtV;
	};

	void LowRank::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		EucRepToExtr(x, egf);
		egf->CopyTo(gf);
		if (prob->GetUseHess())
		{
			Vector *segf = egf->ConstructEmpty();
			egf->CopyTo(segf);
			SharedSpace *Sharedegf = new SharedSpace(segf);
			x->AddToTempData("EGrad", Sharedegf);
		}
		return;
	};

	void LowRank::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		/*compute xix = P_{T_x M} (exix) */
		exix->CopyTo(xix);
		EucRepToExtr(x, xix);

		/*Compute the extra two terms*/
		const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
		Vector *segf = Sharedegf->GetSharedElement();

		LowRankVector *LRetax = dynamic_cast<LowRankVector *> (etax);
		const double *dU = LRetax->GetElement(0)->ObtainReadData();
		const double *dV = LRetax->GetElement(2)->ObtainReadData();

		const SharedSpace *EucRep = segf->ObtainReadTempData("EucRep");
		const double *EucRepptr = EucRep->ObtainReadData();
		const double *M = nullptr;
		blas_sparse_matrix sM;

		/*Obtain the Euclidean representation of the tangent vector,
		if it is sparse, then store in sM, otherwise, store in M.*/
		if (EucRepptr[0]) /*the first entry of EucRepptr is used to indicate if sparse format is used.*/
		{ /*if yes, then use sparse format*/
			sM = static_cast<blas_sparse_matrix> (segf->ObtainReadTempData("SparseMatrix")->ObtainReadData()[0]);
		}
		else /*otherwise, it is a dense matrix.*/
		{
			if (EucRep->Getlength() != 1 + n * m)
			{
				printf("Error: The format of a dense matrix in LowRank::EucHvToHv is not correct!");
			}
			M = EucRepptr + 1;
		}

		/*Compute MV = M * V or sM * V and MTV = M^T * U or sM^T * U */
		double *MdV = new double[m * r], *MTdU = new double[n * r];
		LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
		const double *U = LRx->GetElement(0)->ObtainReadData();
		const double *D = LRx->GetElement(1)->ObtainReadData();
		const double *V = LRx->GetElement(2)->ObtainReadData();

		integer rr = r, mm = m, nn = n;
		if (EucRepptr[0]) /*the first entry of EucRepptr is used to indicate if sparse format is used*/
		{
			for (integer i = 0; i < m * r; i++)
				MdV[i] = 0;
			for (integer i = 0; i < n * r; i++)
				MTdU[i] = 0;
			/*MdV = sM * dV*/
			BLAS_dusmm(blas_colmajor, blas_no_trans, r, 1, sM, dV, n, MdV, m);
			/*MTdU = M^T * dU*/
			BLAS_dusmm(blas_colmajor, blas_trans, r, 1, sM, dU, m, MTdU, n);
		}
		else
		{
			/*MdV = M * dV*/
			dgemm_(GLOBAL::N, GLOBAL::N, &mm, &rr, &nn, &GLOBAL::DONE, const_cast<double *> (M), &mm, const_cast<double *> (dV), &nn, &GLOBAL::DZERO, MdV, &mm);
			/*MTdU = M^T * dU*/
			dgemm_(GLOBAL::T, GLOBAL::N, &nn, &rr, &mm, &GLOBAL::DONE, const_cast<double *> (M), &mm, const_cast<double *> (dU), &mm, &GLOBAL::DZERO, MTdU, &nn);
		}

		double *tmp = new double[r * r + r * ((m > n) ? m : n)];
		
		/*Get LU decomposition of D, which has been done*/
		LUofDinx(x);
		const SharedSpace *SharedSpacetmp = x->ObtainReadTempData("LUofD");
		const double *LUofD = SharedSpacetmp->ObtainReadData();
		integer *P = new integer[r];
		for (integer i = 0; i < r; i++)
			P[i] = static_cast<integer> (LUofD[r*r + i]);

		/*Compute dU_R = (I - U U^T) MdV D^{-1} */
		// U^T MdV
		dgemm_(GLOBAL::T, GLOBAL::N, &rr, &rr, &mm, &GLOBAL::DONE, const_cast<double *> (U), &mm, MdV, &mm, &GLOBAL::DZERO, tmp, &rr);
		// U (U^T MdV)
		dgemm_(GLOBAL::N, GLOBAL::N, &mm, &rr, &rr, &GLOBAL::DONE, const_cast<double *> (U), &mm, tmp, &rr, &GLOBAL::DZERO, tmp + r * r, &mm);
		// (I - U U^T) MdV
		integer length = m * r;
		daxpy_(&length, &GLOBAL::DNONE, tmp + r * r, &GLOBAL::IONE, MdV, &GLOBAL::IONE);

		//X D = (I-U U^T) MdV --> D^T X^T = ((I - U U^T) MdV)^T
		/*Compute the transpose of the right hand side*/
		for (integer i = 0; i < m; i++)
			for (integer j = 0; j < r; j++)
				tmp[r * r + j + i * r] = MdV[i + j * m];
		/*Solve the linear system*/
		integer info;
		dgetrs_(GLOBAL::T, &rr, &mm, const_cast<double *> (LUofD), &rr, P, tmp + r * r, &rr, &info);
		/*Take the transpose*/
		for (integer i = 0; i < m; i++)
			for (integer j = 0; j < r; j++)
				MdV[i + j * m] = tmp[r * r + j + i * r];
		if (info != 0)
			printf("Warning: dgetrs in LowRank::EucHvToHv failed!\n");
		/*dU_xix = dU_xix + dU_R*/
		LowRankVector *LRxix = dynamic_cast<LowRankVector *> (xix);
		LRxix->CopyOnWrite();
		double *dU_xix = LRxix->GetElement(0)->ObtainWritePartialData();
		double *dV_xix = LRxix->GetElement(2)->ObtainWritePartialData();
		daxpy_(&length, &GLOBAL::DONE, MdV, &GLOBAL::IONE, dU_xix, &GLOBAL::IONE);

		/*Compute dU_R = (I - V V^T) MTdU D^{-T} */

		// V^T MTdU
		dgemm_(GLOBAL::T, GLOBAL::N, &rr, &rr, &nn, &GLOBAL::DONE, const_cast<double *> (V), &nn, MTdU, &nn, &GLOBAL::DZERO, tmp, &rr);
		// V (V^T MTdU)
		dgemm_(GLOBAL::N, GLOBAL::N, &nn, &rr, &rr, &GLOBAL::DONE, const_cast<double *> (V), &nn, tmp, &rr, &GLOBAL::DZERO, tmp + r * r, &nn);
		// (I - V V^T) MTdU
		length = n * r;
		daxpy_(&length, &GLOBAL::DNONE, tmp + r * r, &GLOBAL::IONE, MTdU, &GLOBAL::IONE);

		//X D^T = (I - V V^T) MTdU --> D X^T = ((I - V V^T) MTdU)^T

		/*Compute the transpose of the right hand side*/
		for (integer i = 0; i < n; i++)
			for (integer j = 0; j < r; j++)
				tmp[r * r + j + i * r] = MTdU[i + j * n];
		/*Solve the linear system*/
		dgetrs_(GLOBAL::N, &rr, &nn, const_cast<double *> (LUofD), &rr, P, tmp + r * r, &rr, &info);
		/*Take the transpose*/
		for (integer i = 0; i < n; i++)
			for (integer j = 0; j < r; j++)
				MTdU[i + j * n] = tmp[r * r + j + i * r];
		if (info != 0)
			printf("Warning: dgetrs in LowRank::EucHvToHv failed!\n");

		/*dV_xix = dV_xix + dV_R*/
		daxpy_(&length, &GLOBAL::DONE, MTdU, &GLOBAL::IONE, dV_xix, &GLOBAL::IONE);

		delete[] MdV;
		delete[] MTdU;
		delete[] P;
		delete[] tmp;
	};

	blas_sparse_matrix LowRank::ConstructSparseMatrix(Vector *result) const
	{
		const SharedSpace *EucRep = result->ObtainReadTempData("EucRep");
		if (EucRep == nullptr)
			return -1;
		const double *EucRepptr = EucRep->ObtainReadData();
		if (! EucRepptr[0])
			return -1;
		const double *B = nullptr;
		int *ir = nullptr, *jc = nullptr;
		integer nzmax = 0;
		blas_sparse_matrix sM;

		nzmax = static_cast<integer> (EucRepptr[1]);
		if (EucRep->Getlength() != 2 + 3 * nzmax)
		{
			printf("Error: The format of a sparse matrix in LowRank::ConstructSparseMatrix is not correct!");
			return -1;
		}

		B = EucRepptr + 2;
		ir = new int[nzmax];
		jc = new int[nzmax];
		for (integer i = 0; i < nzmax; i++)
		{
			ir[i] = static_cast<int> (EucRepptr[2 + nzmax + i]);
			jc[i] = static_cast<int> (EucRepptr[2 + 2 * nzmax + i]);
		}
		/*Create a sparse matrix using the sparse BLAS library*/
		sM = BLAS_duscr_begin(m, n);
		BLAS_duscr_insert_entries(sM, nzmax, B, ir, jc);
		BLAS_duscr_end(sM);
		if (ir != nullptr)
			delete[] ir;
		if (jc != nullptr)
			delete[] jc;
		return sM;
	};

	void LowRank::EucRepToExtr(Variable *x, Vector *result) const
	{
		const SharedSpace *EucRep = result->ObtainReadTempData("EucRep");
		const double *EucRepptr = EucRep->ObtainReadData();
		const double *M = nullptr;
		blas_sparse_matrix sM;

		/*Obtain the Euclidea representation of the tangent vector,
		if it is sparse, then store in sM, otherwise, store in M.*/
		if (EucRepptr[0]) /*the first entry of EucRepptr is used to indicate if sparse format is used.*/
		{ /*if yes, then use sparse format*/
			sM = ConstructSparseMatrix(result);
		}
		else /*otherwise, it is a dense matrix.*/
		{
			if (EucRep->Getlength() != 1 + n * m)
			{
				printf("Error: The format of a dense matrix in LowRank::EucRepToExtr is not correct!");
			}
			M = EucRepptr + 1;
		}

		/*Compute MV = M * V or sM * V and MTV = M^T * U or sM^T * U */
		double *MV = new double[m * r], *MTU = new double[n * r];
		LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
		const double *U = LRx->GetElement(0)->ObtainReadData();
		const double *D = LRx->GetElement(1)->ObtainReadData();
		const double *V = LRx->GetElement(2)->ObtainReadData();

		integer rr = r, mm = m, nn = n;
		if (EucRepptr[0]) /*the first entry of EucRepptr is used to indicate if sparse format is used*/
		{
			for (integer i = 0; i < m * r; i++)
				MV[i] = 0;
			for (integer i = 0; i < n * r; i++)
				MTU[i] = 0;
			/*MV = sM * V*/
			BLAS_dusmm(blas_colmajor, blas_no_trans, r, 1, sM, V, n, MV, m);
			/*MTU = M^T * U*/
			BLAS_dusmm(blas_colmajor, blas_trans, r, 1, sM, U, m, MTU, n);
		}
		else
		{
			/* MV = M * V */
			dgemm_(GLOBAL::N, GLOBAL::N, &mm, &rr, &nn, &GLOBAL::DONE, const_cast<double *> (M), &mm, const_cast<double *> (V), &nn, &GLOBAL::DZERO, MV, &mm);
			/*MTU = M^T * U*/
			dgemm_(GLOBAL::T, GLOBAL::N, &nn, &rr, &mm, &GLOBAL::DONE, const_cast<double *> (M), &mm, const_cast<double *> (U), &mm, &GLOBAL::DZERO, MTU, &nn);
		}

		LowRankVector *LRresult = dynamic_cast<LowRankVector *> (result);
		LRresult->NewMemoryOnWrite();

		double *dotU = LRresult->GetElement(0)->ObtainWriteEntireData();
		double *dotD = LRresult->GetElement(1)->ObtainWriteEntireData();
		double *dotV = LRresult->GetElement(2)->ObtainWriteEntireData();

		/*compute dotD = U^T MV*/
		dgemm_(GLOBAL::T, GLOBAL::N, &rr, &rr, &mm, &GLOBAL::DONE, const_cast<double *> (U), &mm, MV, &mm, &GLOBAL::DZERO, dotD, &rr);

		integer length = r * m;
		dcopy_(&length, MV, &GLOBAL::IONE, dotU, &GLOBAL::IONE);
		length = r * n;
		dcopy_(&length, MTU, &GLOBAL::IONE, dotV, &GLOBAL::IONE);

		/*compute MV - U dotD*/
		dgemm_(GLOBAL::N, GLOBAL::N, &mm, &rr, &rr, &GLOBAL::DNONE, const_cast<double *> (U), &mm, dotD, &rr, &GLOBAL::DONE, dotU, &mm);

		/*compute MTU - V dotD^T*/
		dgemm_(GLOBAL::N, GLOBAL::T, &nn, &rr, &rr, &GLOBAL::DNONE, const_cast<double *> (V), &nn, dotD, &rr, &GLOBAL::DONE, dotV, &nn);
		
		/*Compute the LU docomposition of D for solving a linear system*/
		LUofDinx(x);
		const SharedSpace *SharedSpacetmp = x->ObtainReadTempData("LUofD");
		const double *tmp = SharedSpacetmp->ObtainReadData();
		integer info;
		integer *P = new integer[r], rsquare = r * r;
		for (integer i = 0; i < r; i++)
			P[i] = static_cast<integer> (tmp[r * r + i]);

		/*Compute dotU = (MV - U dotD) D^{-1} by solving the linear system
			D^T X^T = (MV - U dotD)^T. */
		/*Compute the transpose of the right hand side*/
		double *tmpM = new double[r * ((m > n) ? m : n)];
		for (integer i = 0; i < m; i++)
			for (integer j = 0; j < r; j++)
				tmpM[j + i * r] = dotU[i + j * m];
		/*Solve the linear system*/
		dgetrs_(GLOBAL::T, &rr, &mm, const_cast<double *> (tmp), &rr, P, tmpM, &rr, &info);
		/*dotU is the transpose of the solution*/
		for (integer i = 0; i < m; i++)
			for (integer j = 0; j < r; j++)
				dotU[i + j * m] = tmpM[j + i * r];

		if (info != 0)
			printf("Warning: dgetrs in LowRank::EucRepToExtr failed!\n");

		/*Compute dotV = (MTU - V dotD) D^{-T} by solving the linear system
			D X^T = (MTU - V dotD)^T */

		/*Compute the transpose of the right hand side*/
		for (integer i = 0; i < n; i++)
			for (integer j = 0; j < r; j++)
				tmpM[j + i * r] = dotV[i + j * n];
		/*Solve the linear system*/
		dgetrs_(GLOBAL::N, &rr, &nn, const_cast<double *> (tmp), &rr, P, tmpM, &rr, &info);
		/*dotV is the transpose of the solution*/
		for (integer i = 0; i < n; i++)
			for (integer j = 0; j < r; j++)
				dotV[i + j * n] = tmpM[j + i * r];

		if (info != 0)
			printf("Warning: dgetrs in LowRank::ProjectionM failed!\n");

		delete[] tmpM;

		if (EucRepptr[0]) // if it is sparse, then attached the sparse matrix to result.
		{
			SharedSpace *SparseMatrix = new SharedSpace(1, 1);
			SparseMatrix->ObtainWriteEntireData()[0] = static_cast<blas_sparse_matrix> (sM);
			result->AddToTempData("SparseMatrix", SparseMatrix);
		}

		delete[] MV;
		delete[] MTU;
		delete[] P;
	};

	void LowRank::LUofDinx(Variable *x) const
	{
		if (!x->TempDataExist("LUofD"))
		{
			LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
			const double *D = LRx->GetElement(1)->ObtainReadData();

			/*Compute the LU docomposition of D for solving a linear system*/
			SharedSpace *SharedSpacetmp = new SharedSpace(1, r*r + r);
			double *tmp = SharedSpacetmp->ObtainWriteEntireData();
			integer info;
			integer *P = new integer[r], rsquare = r * r, rr = r;
			dcopy_(&rsquare, const_cast<double *> (D), &GLOBAL::IONE, tmp, &GLOBAL::IONE);
			// LU decomposion for D, D = P * L * U, L and U are stored in D, the permutation matrix is in P
			// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
			dgetrf_(&rr, &rr, tmp, &rr, P, &info);
			for (integer i = 0; i < r; i++)
				tmp[r * r + i] = static_cast<double> (P[i]);
			if (info != 0)
				printf("Warning: dgetrs in LowRank::LUofDinx failed!\n");
			delete[] P;
			x->AddToTempData("LUofD", SharedSpacetmp);
		}
	};

	void LowRank::ExtrToEucRep(Variable *x, Vector *result) const
	{
		LowRankVariable *LRx = dynamic_cast<LowRankVariable *> (x);
		const double *U = LRx->GetElement(0)->ObtainReadData();
		const double *D = LRx->GetElement(1)->ObtainReadData();
		const double *V = LRx->GetElement(2)->ObtainReadData();
		LowRankVector *LRresult = dynamic_cast<LowRankVector *> (result);
		const double *dotU = LRresult->GetElement(0)->ObtainReadData();
		const double *dotD = LRresult->GetElement(1)->ObtainReadData();
		const double *dotV = LRresult->GetElement(2)->ObtainReadData();

		SharedSpace *EucRep = new SharedSpace(1, 1 + m * n);
		double *EucRepptr = EucRep->ObtainWriteEntireData();
		EucRepptr[0] = 0; //use dense matrix format

		double *tmp = new double[r * n];

		integer rr = r, nn = n, mm = m;
		/*tmp = dotD * V^T */
		dgemm_(GLOBAL::N, GLOBAL::T, &rr, &nn, &rr, &GLOBAL::DONE, const_cast<double *> (dotD), &rr, const_cast<double *> (V), &nn, &GLOBAL::DZERO, tmp, &rr);
		/*result = U * dotD * V^T */
		dgemm_(GLOBAL::N, GLOBAL::N, &mm, &nn, &rr, &GLOBAL::DONE, const_cast<double *> (U), &mm, tmp, &rr, &GLOBAL::DZERO, EucRepptr + 1, &mm);

		/*tmp = D * V^T */
		dgemm_(GLOBAL::N, GLOBAL::T, &rr, &nn, &rr, &GLOBAL::DONE, const_cast<double *> (D), &rr, const_cast<double *> (V), &nn, &GLOBAL::DZERO, tmp, &rr);
		/*result = (U * dotD * V^T) + dotU * D * V^T */
		dgemm_(GLOBAL::N, GLOBAL::N, &mm, &nn, &rr, &GLOBAL::DONE, const_cast<double *> (dotU), &mm, tmp, &rr, &GLOBAL::DONE, EucRepptr + 1, &mm);

		/*tmp = D * dotV^T */
		dgemm_(GLOBAL::N, GLOBAL::T, &rr, &nn, &rr, &GLOBAL::DONE, const_cast<double *> (D), &rr, const_cast<double *> (dotV), &nn, &GLOBAL::DZERO, tmp, &rr);
		/*result = (U * dotD * V^T + dotU * D * V^T) + U * D * dotV^T */
		dgemm_(GLOBAL::N, GLOBAL::N, &mm, &nn, &rr, &GLOBAL::DONE, const_cast<double *> (U), &mm, tmp, &rr, &GLOBAL::DONE, EucRepptr + 1, &mm);

		delete[] tmp;
		result->AddToTempData("EucRep", EucRep);
	}

}; /*end of ROPTLIB namespace*/
