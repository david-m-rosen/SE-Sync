
#include "Manifolds/SPDManifold/SPDManifold.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDManifold::SPDManifold(integer inn)
	{
		n = inn;
		IsIntrApproach = true;
		HasHHR = false;
		UpdBetaAlone = false;
		HasLockCon = false;
		name.assign("SPDManifold");
		IntrinsicDim = n * (n + 1) / 2;
		ExtrinsicDim = n * n;
		EMPTYEXTR = new SPDVector(n, n);
		EMPTYINTR = new SPDVector(IntrinsicDim, 1);
	};

	SPDManifold::~SPDManifold(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	void SPDManifold::CheckParams(void) const
	{
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		printf("row           :%15d,\t", n);
		printf("col           :%15d\n", n);
	};

	void SPDManifold::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		const double *xptr = x->ObtainReadData();
		const double *egfptr = egf->ObtainReadData();
		double *gfptr = gf->ObtainWriteEntireData();
		integer dim = n;
		double *tmp = new double[dim * dim];
		dgemm_(GLOBAL::T, GLOBAL::N, &dim, &dim, &dim, &GLOBAL::DONE, const_cast<double *> (xptr), &dim, const_cast<double *> (egfptr), &dim,
			&GLOBAL::DZERO, tmp, &dim);
		dgemm_(GLOBAL::N, GLOBAL::N, &dim, &dim, &dim, &GLOBAL::DONE, const_cast<double *> (tmp), &dim, const_cast<double *> (xptr), &dim,
			&GLOBAL::DZERO, gfptr, &dim);
		delete[] tmp;
	};

	void SPDManifold::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const
	{
		printf("warning:SPDManifold::EucHvToHv has not been done!\n");
		exix->CopyTo(xix);
	};

	void SPDManifold::CholeskyRepresentation(Variable *x) const
	{
		const double *xM = x->ObtainReadData();
		Variable *L = x->ConstructEmpty();
		SharedSpace *SharedL = new SharedSpace(L);
		double *LM = L->ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i; j < n; j++)
			{
				LM[i + j * n] = 0;
				LM[j + i * n] = xM[j + i * n];
			}
		}

		integer info, N = n;
		dpotrf_(GLOBAL::L, &N, LM, &N, &info);
		x->AddToTempData("L", SharedL);
		if (info != 0)
		{
			printf("Warning: SPDManifold::CholeskyRepresentation fails with info:%d!\n", info);
		}
	};

	void SPDManifold::ExtrProjection(Variable *x, Vector *etax, Vector *result) const
	{
		const double *etaxTV = etax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
		{
			resultTV[i + i * n] = etaxTV[i + i * n];
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[i + j * n] = (etaxTV[i + j * n] + etaxTV[j + i * n]) * 0.5;
				resultTV[j + i * n] = resultTV[i + j * n];
			}
		}
	};

	void SPDManifold::ObtainIntr(Variable *x, Vector *etax, Vector *result) const
	{//L^{-1} etax L^{-T}, where x = L L^T, etax is assumed to be a symmetric matrix.
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();
		double *E = new double[n * n];
		integer length = n * n, N = n, info;
		dcopy_(&length, const_cast<double *> (etax->ObtainReadData()), &GLOBAL::IONE, E, &GLOBAL::IONE);
		/*Solve the linear system L X = E, i.e., X = L^{-1} E. The solution X is stored in E
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html*/
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, E, &N, &info);
		if (info != 0)
		{
			printf("warning: SPDManifold::ObtainIntr fails with info:%d!\n", info);
		}

		/*E <-- E^T*/
		double tmp;
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				tmp = E[i + j * n];
				E[i + j * n] = E[j + i * n];
				E[j + i * n] = tmp;
			}
		}
		/*Solve the linear system L X = E, i.e., X = L^{-1} E. The solution X is stored in E
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, E, &N, &info);
		if (info != 0)
		{
			printf("warning: SPDManifold::ObtainIntr fails with info:%d!\n", info);
		}
		/*We don't have to do: E <-- E^T, since E is symmetric*/
		double *resultTV = result->ObtainWriteEntireData();
		integer idx = 0;
		double r2 = sqrt(2.0);
		for (integer i = 0; i < n; i++)
		{
			resultTV[idx] = E[i + i * n];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[idx] = E[j + i * n] * r2;
				idx++;
			}
		}

		delete[] E;
	};

	void SPDManifold::ObtainExtr(Variable *x, Vector *intretax, Vector *result) const
	{
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();
		const double *intretaxTV = intretax->ObtainReadData();
		double *resultTV = result->ObtainWriteEntireData();

		integer idx = 0;
		double r2 = sqrt(2.0);
		for (integer i = 0; i < n; i++)
		{
			resultTV[i + i * n] = intretaxTV[idx];
			idx++;
		}

		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[j + i * n] = intretaxTV[idx] / r2;
				resultTV[i + j * n] = resultTV[j + i * n];
				idx++;
			}
		}

		double *E = new double[n * n];
		integer N = n;

		/*E <-- L resultTV */
		dgemm_(GLOBAL::N, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, const_cast<double *> (L), &N, resultTV, &N,
			&GLOBAL::DZERO, E, &N);
		/*resultTV <-- E L^T */
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, E, &N, const_cast<double *> (L), &N,
			&GLOBAL::DZERO, resultTV, &N);

		delete[] E;
	};

	double SPDManifold::Dist(Variable *x1, Variable *x2) const
	{
		// dis(x, y) = |log(Lx^{-1} y Lx^{-T})|_F
		if (!x1->TempDataExist("L"))
		{
			CholeskyRepresentation(x1);
		}
		const SharedSpace *SharedLx1 = x1->ObtainReadTempData("L");
		Variable *LElementx1 = SharedLx1->GetSharedElement();
		const double *Lx1 = LElementx1->ObtainReadData();

		if (!x2->TempDataExist("L"))
		{
			CholeskyRepresentation(x2);
		}
		const SharedSpace *SharedLx2 = x2->ObtainReadTempData("L");
		Variable *LElementx2 = SharedLx2->GetSharedElement();
		const double *Lx2 = LElementx2->ObtainReadData();

		integer N = n, info, length = n * n;
		double *LxLy = new double[N * N];

		dcopy_(&length, const_cast<double *> (Lx2), &GLOBAL::IONE, LxLy, &GLOBAL::IONE);
		/*Solve the linear system Lx X = Ly, i.e., X = Lx^{-1}Ly. The solution X is stored in LxLy*/
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (Lx1), &N, LxLy, &N, &info);

		/* compute temp = LxLy*(LxLy)^T  */
		double *temp = new double[N * N];
		dgemm_(GLOBAL::N, GLOBAL::T, &N, &N, &N, &GLOBAL::DONE, LxLy, &N, LxLy, &N, &GLOBAL::DZERO, temp, &N);

		double *eigenvalues = new double[n + n * n];
		double *eigenvectors = eigenvalues + n;
		Matrix E(eigenvalues, n, 1), V(eigenvectors, n, n);
		Matrix MMt(temp, n, n);
		Matrix::EigenSymmetricM(GLOBAL::L, MMt, E, V);

		double result = 0.0;
		for (integer i = 0; i < N; i++)
		{
			result = result + log(eigenvalues[i])*log(eigenvalues[i]);
		}

		result = sqrt(result);

		delete[] eigenvalues;
		delete[] temp;
		delete[] LxLy;
		return result;
	};

	void SPDManifold::Retraction(Variable *x, Vector *etax, Variable *result) const
	{
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();

		/*Compute the extrinsic representation*/
		Vector *exetax = nullptr;
		const double *etaxTV = nullptr;
		if (IsIntrApproach)
		{
			exetax = EMPTYEXTR->ConstructEmpty();
			ObtainExtr(x, etax, exetax);
			etaxTV = exetax->ObtainReadData();
		}
		else
		{
			etaxTV = etax->ObtainReadData();
		}

		integer N = n, info;
		double *LiE = new double[N * N];
		integer length = n * n;
		dcopy_(&length, const_cast<double *> (etaxTV), &GLOBAL::IONE, LiE, &GLOBAL::IONE);
		/*Solve the linear system L X = E, i.e., X = L^{-1} E. The solution X is stored in LiE
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, LiE, &N, &info);
		if (info != 0)
		{
			printf("warning: SPDManifold::Retraction fails with info:%d!\n", info);
		}
		double *resultTV = result->ObtainWriteEntireData();
		/*Compute result = LiE^T LiE = E L^{-T} L^{-1} E*/
		dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, LiE, &N, LiE, &N, &GLOBAL::DZERO, resultTV, &N);
		delete[] LiE;
		double scalar = 0.5;
		/*result <-- 0.5 * result*/
		dscal_(&length, &scalar, resultTV, &GLOBAL::IONE);

		/*result <-- etax + result*/
		daxpy_(&length, &GLOBAL::DONE, const_cast<double *> (etaxTV), &GLOBAL::IONE, resultTV, &GLOBAL::IONE);

		const double *xM = x->ObtainReadData();
		/*result <-- x + result*/
		daxpy_(&length, &GLOBAL::DONE, const_cast<double *> (xM), &GLOBAL::IONE, resultTV, &GLOBAL::IONE);
		delete exetax;

		if (!result->TempDataExist("L"))
		{
			CholeskyRepresentation(result);
		}
	};

	void SPDManifold::DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const
	{ // xix + 0.5 xix x^{-1} etax + 0.5 etax x^{-1} xix
		if (!x->TempDataExist("L"))
		{
			CholeskyRepresentation(x);
		}
		const SharedSpace *SharedL = x->ObtainReadTempData("L");
		Variable *LElement = SharedL->GetSharedElement();
		const double *L = LElement->ObtainReadData();

		/*Compute the extrinsic representation*/
		Vector *exetax = EMPTYEXTR->ConstructEmpty();
		Vector *exxix = EMPTYEXTR->ConstructEmpty();
		ObtainExtr(x, etax, exetax);
		ObtainExtr(x, xix, exxix);
		double *LiE = new double[2 * n * n];
		double *LiX = LiE + n * n;
		const double *exetaxTV = exetax->ObtainReadData();
		const double *exxixTV = exxix->ObtainReadData();
		integer length = n * n, N = n, info;
		dcopy_(&length, const_cast<double*> (exetaxTV), &GLOBAL::IONE, LiE, &GLOBAL::IONE);
		dcopy_(&length, const_cast<double*> (exxixTV), &GLOBAL::IONE, LiX, &GLOBAL::IONE);
		delete exetax;

		/*Solve the linear system L X = Eta, i.e., X = L^{-1} Eta. The solution X is stored in LiE
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, LiE, &N, &info);

		/*Solve the linear system L X = Xix, i.e., X = L^{-1} Xix. The solution X is stored in LiX
		Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
		dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &N, &N, const_cast<double *> (L), &N, LiX, &N, &info);

		Vector *exresult = EMPTYEXTR->ConstructEmpty();
		double *resultTV = exresult->ObtainWriteEntireData();
		/*resultTV <-- etax L^{-T} L^{-1} xix*/
		dgemm_(GLOBAL::T, GLOBAL::N, &N, &N, &N, &GLOBAL::DONE, LiE, &N, LiX, &N, &GLOBAL::DZERO, resultTV, &N);
		for (integer i = 0; i < n; i++)
		{
			for (integer j = i + 1; j < n; j++)
			{
				resultTV[j + i * n] = (resultTV[j + i * n] + resultTV[i + j * n]) * 0.5;
				resultTV[i + j * n] = resultTV[j + i * n];
			}
		}
		delete[] LiE;

		daxpy_(&length, &GLOBAL::DONE, const_cast<double*> (exxixTV), &GLOBAL::IONE, resultTV, &GLOBAL::IONE);
		delete exxix;
		ObtainIntr(y, exresult, result);
		delete exresult;

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

	void SPDManifold::coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const
	{
		printf("SPDManifold::coTangentVector has not been done!\n");
		xiy->CopyTo(result);
	};

	double SPDManifold::Beta(Variable *x, Vector *etax) const
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
}; /*end of ROPTLIB namespace*/
