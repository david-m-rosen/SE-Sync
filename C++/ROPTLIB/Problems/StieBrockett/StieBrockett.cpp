
#include "Problems/StieBrockett/StieBrockett.h"

/*Define the namespace*/
namespace ROPTLIB{

	StieBrockett::StieBrockett(double *inB, double *inD, integer inn, integer inp)
	{
		B = inB;
		D = inD;
		n = inn;
		p = inp;
	};

	StieBrockett::~StieBrockett(void)
	{
	};

	double StieBrockett::f(Variable *x) const
	{
		const double *xxM = x->ObtainReadData();
		Vector *BxD = x->ConstructEmpty();
		SharedSpace *Temp = new SharedSpace(BxD);
		double *temp = BxD->ObtainWriteEntireData();
		double result = 0;

		Matrix MB(B, n, n), MxxM(xxM, n, p), Mtemp(temp, n, p);
		// temp = B * xxM, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		Matrix::DGEMM(1, MB, false, MxxM, false, 0, Mtemp);

		for (integer i = 0; i < p; i++)
		{
			// temp(i * n : i * n + N - 1) <- D[i] * temp(i * n : i * n + N - 1),
			// details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(const_cast<integer *> (&n), &D[i], temp + i * n, &GLOBAL::IONE);
		}
		integer length = n * p;
		// output temp(:)^T * xxM(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		result = ddot_(&length, temp, &GLOBAL::IONE, const_cast<double *> (xxM), &GLOBAL::IONE);
		if (UseGrad)
		{
			x->AddToTempData("BxD", Temp);
		}
		else
		{
			delete Temp;
		}
		return result;
	};

	//void StieBrockett::RieGrad(Variable *x, Vector *gf) const
	//{
	//	const double *xxM = x->ObtainReadData();
	//	const SharedSpace *SharedTemp = x->ObtainReadTempData("BxD");
	//	Vector *Temp = SharedTemp->GetSharedElement();
	//	const double *xBxD = Temp->ObtainReadData();
	//	double *gfTV = gf->ObtainWriteEntireData();
	//	SharedSpace *symxtBxN = new SharedSpace(2, p, p);
	//	double *symxtBxNptr = symxtBxN->ObtainWriteEntireData();
	//
	//	char *transn = "n", *transt = "t";
	//	double one = 1, zero = 0;
	//	integer inc = 1, N = n, P = p, Length = N * P;
	//	dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (xxM), &N, const_cast<double *> (xBxD), &N, &zero, symxtBxNptr, &P);
	//	for (integer i = 0; i < p; i++)
	//	{
	//		for (integer j = i + 1; j < p; j++)
	//		{
	//			symxtBxNptr[i + j * p] += symxtBxNptr[j + i * p];
	//			symxtBxNptr[i + j * p] /= 2.0;
	//			symxtBxNptr[j + i * p] = symxtBxNptr[i + j * p];
	//		}
	//	}
	//	dcopy_(&Length, const_cast<double *> (xBxD), &inc, gfTV, &inc);
	//	double negone = -1;
	//	dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (xxM), &N, symxtBxNptr, &P, &one, gfTV, &N);
	//	double two = 2;
	//	dscal_(&Length, &two, gfTV, &inc);
	//	if (UseHess)
	//	{
	//		x->AddToTempData("symxtBxD", symxtBxN);
	//	}
	//	else
	//	{
	//		delete symxtBxN;
	//	}
	//};
	//
	//void StieBrockett::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	//{
	//	if (!x->TempDataExist("symxtBxD"))
	//	{
	//		const double *xxM = x->ObtainReadData();
	//		const SharedSpace *Temp = x->ObtainReadTempData("BxD");
	//		const double *xBxD = Temp->ObtainReadData();
	//		SharedSpace *symxtBxN = new SharedSpace(2, p, p);
	//		double *symxtBxNptr = symxtBxN->ObtainWriteEntireData();
	//
	//		char *transn = "n", *transt = "t";
	//		double one = 1, zero = 0;
	//		integer inc = 1, N = n, P = p, Length = N * P;
	//		dgemm_(transt, transn, &P, &P, &N, &one, const_cast<double *> (xxM), &N, const_cast<double *> (xBxD), &N, &zero, symxtBxNptr, &P);
	//		for (integer i = 0; i < p; i++)
	//		{
	//			for (integer j = i + 1; j < p; j++)
	//			{
	//				symxtBxNptr[i + j * p] += symxtBxNptr[j + i * p];
	//				symxtBxNptr[i + j * p] /= 2.0;
	//				symxtBxNptr[j + i * p] = symxtBxNptr[i + j * p];
	//			}
	//		}
	//		x->AddToTempData("symxtBxD", symxtBxN);
	//	}
	//
	//	const SharedSpace *symxtBxD = x->ObtainReadTempData("symxtBxD");
	//	const double *DsymxtBxD = symxtBxD->ObtainReadData();
	//	const double *xM = x->ObtainReadData();
	//	const double *etaxTV = etax->ObtainReadData();
	//	double *xixTV = xix->ObtainWriteEntireData();
	//
	//	if (xixTV == etaxTV)
	//	{
	//		printf("Error in RieHessianEta!\n");
	//		exit(0);
	//	}
	//	char *transn = "n";
	//	integer N = n, P = p, inc = 1, Length = N * P;
	//	double one = 1, zero = 0, negone = -1, two = 2;
	//	dgemm_(transn, transn, &N, &P, &N, &one, B, &N, const_cast<double *> (etaxTV), &N, &zero, xixTV, &N);
	//	for (integer i = 0; i < p; i++)
	//	{
	//		dscal_(&N, &D[i], xixTV + i * n, &inc);
	//	}
	//
	//	dgemm_(transn, transn, &N, &P, &P, &negone, const_cast<double *> (etaxTV), &N, const_cast<double *> (DsymxtBxD), &P, &one, xixTV, &N);
	//	dscal_(&Length, &two, xixTV, &inc);
	//	Stiefel *StieDomain = dynamic_cast<Stiefel *> (Domain);
	//	StieDomain->ExtrProjection(x, xix, xix);
	//};

	void StieBrockett::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *Temp = x->ObtainReadTempData("BxD");
		Vector *BxD = Temp->GetSharedElement();
		Domain->ScaleTimesVector(x, 2.0, BxD, egf);
	};

	void StieBrockett::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const double *etaxTV = etax->ObtainReadData();
		double *exixTV = exix->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n");
		integer N = n, P = p, inc = 1, Length = N * P;
		double one = 1, zero = 0, two = 2;
		// exxiTV <- B * etaxTV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &N, &one, B, &N, const_cast<double *> (etaxTV), &N, &zero, exixTV, &N);
		for (integer i = 0; i < p; i++)
		{
			// exixTV(i * n : i * n + N - 1) <- D[i] * exixTV(i * n : i * n + N - 1),
			// details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &D[i], exixTV + i * n, &inc);
		}
		Domain->ScaleTimesVector(x, 2.0, exix, exix);
	};
}; /*end of ROPTLIB namespace*/
