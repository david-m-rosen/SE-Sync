
#include "Problems/StieSumBrockett/StieSumBrockett.h"

/*Define the namespace*/
namespace ROPTLIB{

	StieSumBrockett::StieSumBrockett(double *inB1, double *inD1, double *inB2, double *inD2, double *inB3,
		double *inD3, integer inn, integer inp, integer inm, integer inq)
	{
		B1 = inB1;
		D1 = inD1;
		B2 = inB2;
		D2 = inD2;
		B3 = inB3;
		D3 = inD3;
		n = inn;
		p = inp;
		m = inm;
		q = inq;
	};

	StieSumBrockett::~StieSumBrockett(void)
	{
	};

	double StieSumBrockett::f(Variable *x) const
	{
		const double *xX1 = x->ObtainReadData();
		const double *xX2 = xX1 + n * p;
		const double *xX3 = xX2 + n * p;

		ProductElement *prodx = dynamic_cast<ProductElement *> (x);
		Vector *BxD1 = prodx->GetElement(0)->ConstructEmpty();
		SharedSpace *Temp1 = new SharedSpace(BxD1);
		double *temp1 = BxD1->ObtainWriteEntireData();
		double result = 0;

		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		integer inc = 1, N = n, P = p;
		// temp1 <- B1 * xX1, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &N, &one, B1, &N, const_cast<double *> (xX1), &N, &zero, temp1, &N);
		for (integer i = 0; i < p; i++)
		{
			// temp(:, i) <- D1[i] * temp(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &D1[i], temp1 + i * n, &inc);
		}
		integer length = N * P;
		// compute temp1(:)^T xX1(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		result += ddot_(&length, temp1, &inc, const_cast<double *> (xX1), &inc);
		if (UseGrad)
		{
			x->AddToTempData("BxD1", Temp1);
		}
		else
		{
			delete Temp1;
		}

		Vector *BxD2 = prodx->GetElement(1)->ConstructEmpty();
		SharedSpace *Temp2 = new SharedSpace(BxD2);
		double *temp2 = BxD2->ObtainWriteEntireData();

		// temp2 <- B2 * xX2, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &N, &one, B2, &N, const_cast<double *> (xX2), &N, &zero, temp2, &N);
		for (integer i = 0; i < p; i++)
		{
			// temp2(:, i) <- D2[i] * temp2(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &D2[i], temp2 + i * n, &inc);
		}
		// compute temp2(:)^T xX2(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		result += ddot_(&length, temp2, &inc, const_cast<double *> (xX2), &inc);
		if (UseGrad)
		{
			x->AddToTempData("BxD2", Temp2);
		}
		else
		{
			delete Temp2;
		}

		Vector *BxD3 = prodx->GetElement(2)->ConstructEmpty();
		SharedSpace *Temp3 = new SharedSpace(BxD3);
		double *temp3 = BxD3->ObtainWriteEntireData();
		integer M = m, Q = q;
		length = M * Q;
		// temp2 <- B2 * xX2, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &M, &Q, &M, &one, B3, &M, const_cast<double *> (xX3), &M, &zero, temp3, &M);
		for (integer i = 0; i < q; i++)
		{
			// temp3(:, i) <- D3[i] * temp3(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&M, &D3[i], temp3 + i * m, &inc);
		}
		// compute temp3(:)^T xX3(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		result += ddot_(&length, temp3, &inc, const_cast<double *> (xX3), &inc);
		if (UseGrad)
		{
			x->AddToTempData("BxD3", Temp3);
		}
		else
		{
			delete Temp3;
		}

		return result;
	};

	void StieSumBrockett::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *Temp1 = x->ObtainReadTempData("BxD1");
		const SharedSpace *Temp2 = x->ObtainReadTempData("BxD2");
		const SharedSpace *Temp3 = x->ObtainReadTempData("BxD3");
		Vector *BxD1 = Temp1->GetSharedElement();
		Vector *BxD2 = Temp2->GetSharedElement();
		Vector *BxD3 = Temp3->GetSharedElement();

		ProductElement *prodegf = dynamic_cast<ProductElement *> (egf);
		ProductElement *prodx = dynamic_cast<ProductElement *> (x);
		prodegf->NewMemoryOnWrite();

		ProductManifold *ProdDomain = dynamic_cast<ProductManifold *> (Domain);

		ProdDomain->GetManifold(0)->ScaleTimesVector(prodx->GetElement(0), 2.0, BxD1, prodegf->GetElement(0));
		ProdDomain->GetManifold(0)->ScaleTimesVector(prodx->GetElement(1), 2.0, BxD2, prodegf->GetElement(1));
		ProdDomain->GetManifold(1)->ScaleTimesVector(prodx->GetElement(2), 2.0, BxD3, prodegf->GetElement(2));
	};

	void StieSumBrockett::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		ProductElement *prodx = dynamic_cast<ProductElement *> (x);
		ProductElement *prodetax = dynamic_cast<ProductElement *> (etax);
		ProductElement *prodexix = dynamic_cast<ProductElement *> (exix);
		prodexix->NewMemoryOnWrite();
		ProductManifold *ProdDomain = dynamic_cast<ProductManifold *> (Domain);

		const double *etax1TV = prodetax->GetElement(0)->ObtainReadData();
		double *exix1TV = prodexix->GetElement(0)->ObtainWriteEntireData();
		char *transn = const_cast<char *> ("n");
		integer N = n, P = p, inc = 1, Length = N * P;
		double one = 1, zero = 0, negone = -1, two = 2;
		// exix1TV <- B1 * etax1TV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &N, &one, B1, &N, const_cast<double *> (etax1TV), &N, &zero, exix1TV, &N);
		for (integer i = 0; i < p; i++)
		{
			// exix1TV(:, i) <- D1[i] * exix1TV(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &D1[i], exix1TV + i * n, &inc);
		}
		ProdDomain->GetManifold(0)->ScaleTimesVector(prodx->GetElement(0), 2.0, prodexix->GetElement(0), prodexix->GetElement(0));

		const double *etax2TV = prodetax->GetElement(1)->ObtainReadData();
		double *exix2TV = prodexix->GetElement(1)->ObtainWriteEntireData();
		// exix2TV <- B2 * etax2TV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &N, &one, B2, &N, const_cast<double *> (etax2TV), &N, &zero, exix2TV, &N);
		for (integer i = 0; i < p; i++)
		{
			// exix2TV(:, i) <- D2[i] * exix2TV(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &D2[i], exix2TV + i * n, &inc);
		}
		ProdDomain->GetManifold(0)->ScaleTimesVector(prodx->GetElement(1), 2.0, prodexix->GetElement(1), prodexix->GetElement(1));

		const double *etax3TV = prodetax->GetElement(2)->ObtainReadData();
		double *exix3TV = prodexix->GetElement(2)->ObtainWriteEntireData();
		integer M = m, Q = q;
		Length = N * P;
		// exix3TV <- B3 * etax3TV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &M, &Q, &M, &one, B3, &M, const_cast<double *> (etax3TV), &M, &zero, exix3TV, &M);
		for (integer i = 0; i < q; i++)
		{
			// exix3TV(:, i) <- D3[i] * exix3TV(:, i), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&M, &D3[i], exix3TV + i * m, &inc);
		}
		ProdDomain->GetManifold(1)->ScaleTimesVector(prodx->GetElement(2), 2.0, prodexix->GetElement(2), prodexix->GetElement(2));
	};
}; /*end of ROPTLIB namespace*/
