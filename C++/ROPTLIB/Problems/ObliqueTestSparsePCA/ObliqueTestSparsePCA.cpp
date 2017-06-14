
#include "Problems/ObliqueTestSparsePCA/ObliqueTestSparsePCA.h"

/*Define the namespace*/
namespace ROPTLIB{

	ObliqueTestSparsePCA::ObliqueTestSparsePCA(double *inB, double *inDsq, double inmu, double ineps, integer inp, integer inn, integer inr)
	{
		B = inB;
		Dsq = inDsq;
		p = inp;
		n = inn;
		r = inr;
		mu = inmu;
		epsilon = ineps;
	};

	ObliqueTestSparsePCA::~ObliqueTestSparsePCA(void)
	{
	};

	double ObliqueTestSparsePCA::f(Variable *x) const
	{
		const double *xptr = x->ObtainReadData();
		double epsilonsq = epsilon * epsilon;

		double result = 0;
		for (integer i = 0; i < p * r; i++)
		{
			result += sqrt(xptr[i] * xptr[i] + epsilonsq) - epsilon;
		}
		double *BtX = new double[n * r];

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		integer N = n, P = p, R = r;
		double one = 1.0, zero = 0.0;
		// BtX <- B^T * xptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &N, &R, &P, &one, B, &P, const_cast<double *> (xptr), &P, &zero, BtX, &N);
		SharedSpace *SharedBBtX = new SharedSpace(2, p, r);
		double *BBtX = SharedBBtX->ObtainWriteEntireData();
		// BBtX = B * BtX, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &R, &N, &one, B, &P, BtX, &N, &zero, BBtX, &P);

		SharedSpace *SharedXtBBtX = new SharedSpace(2, r, r);
		double *XtBBtX = SharedXtBBtX->ObtainWriteEntireData();
		// XtBBtX = xptr^T * BBtX, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &R, &R, &P, &one, const_cast<double *> (xptr), &P, BBtX, &P, &zero, XtBBtX, &R);

		for (integer i = 0; i < r; i++)
		{
			result += mu * (XtBBtX[i + i * r] - Dsq[i]) * (XtBBtX[i + i * r] - Dsq[i]);
			for (integer j = i + 1; j < r; j++)
			{
				result += mu * 2 * XtBBtX[i + j * r] * XtBBtX[i + j * r];
			}
		}
		//integer rsqr = r * r, inc = 1;
		//result += mu * ddot_(&rsqr, XtBBtX, &inc, Dsq, &inc);

		x->AddToTempData("BBtX", SharedBBtX);
		x->AddToTempData("XtBBtX", SharedXtBBtX);
		delete[] BtX;
		return result;
	};

	void ObliqueTestSparsePCA::EucGrad(Variable *x, Vector *egf) const
	{
		const double *xptr = x->ObtainReadData();
		double *egfptr = egf->ObtainWriteEntireData();
		double epsilonsq = epsilon * epsilon;

		for (integer i = 0; i < p * r; i++)
		{
			egfptr[i] = xptr[i] / sqrt(xptr[i] * xptr[i] + epsilonsq);
		}

		const SharedSpace *SharedXtBBtX = x->ObtainReadTempData("XtBBtX");
		const double *XtBBtX = SharedXtBBtX->ObtainReadData();

		SharedSpace *SharedXtBBtXmDsq = new SharedSpace(2, r, r);
		double *XtBBtXmDsq = SharedXtBBtXmDsq->ObtainWriteEntireData();
		integer rsq = r * r, inc = 1;
		// XtBBtXmDsq <- XtBBtX, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&rsq, const_cast<double *> (XtBBtX), &inc, XtBBtXmDsq, &inc);
		for (integer i = 0; i < r; i++)
		{
			XtBBtXmDsq[i + i * r] -= Dsq[i];
		}

		const SharedSpace *SharedBBtX = x->ObtainReadTempData("BBtX");
		const double *BBtX = SharedBBtX->ObtainReadData();
		integer P = p, R = r;
		char *transn = const_cast<char *> ("n");
		double one = 1.0, fourmu = 4 * mu;
		// egfptr = egfptr + 4 * mu * BBtX * XtBBtXmDsq, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &R, &R, &fourmu, const_cast<double *> (BBtX), &P, XtBBtXmDsq, &R, &one, egfptr, &P);
		x->AddToTempData("XtBBtXmDsq", SharedXtBBtXmDsq);
	};

	void ObliqueTestSparsePCA::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		//etax->CopyTo(exix);
		double epsilonsq = epsilon * epsilon;
		const double *etaxTV = etax->ObtainReadData();
		const double *xptr = x->ObtainReadData();
		double *exixTV = exix->ObtainWriteEntireData();
		for (integer i = 0; i < p * r; i++)
		{
			exixTV[i] = etaxTV[i] * epsilonsq / pow(xptr[i] * xptr[i] + epsilonsq, 1.5);
		}

		const SharedSpace *SharedXtBBtXmDsq = x->ObtainReadTempData("XtBBtXmDsq");
		const double *XtBBtXmDsq = SharedXtBBtXmDsq->ObtainReadData();

		double *etaXtBBtXmDsq = new double[p * r + n * r];
		double *BtetaXtBBtXmDsq = etaXtBBtXmDsq + p * r;

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		integer N = n, P = p, R = r, inc = 1;
		double one = 1.0, zero = 0.0;
		// etaXtBBtXmDsq = etaxTV * XtBBtXmDsq, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &R, &R, &one, const_cast<double *> (etaxTV), &P, const_cast<double *> (XtBBtXmDsq), &R, &zero, etaXtBBtXmDsq, &P);
		// BtetaXtBBtXmDsq = B^T * etaXtBBtXmDsq, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &N, &R, &P, &one, B, &P, etaXtBBtXmDsq, &P, &zero, BtetaXtBBtXmDsq, &N);
		double fourmu = 4 * mu;
		// exixTV = exixTV + 4 * mu * B * BtetaXtBBtXmDsq, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &R, &N, &fourmu, B, &P, BtetaXtBBtXmDsq, &N, &one, exixTV, &P);
		delete[] etaXtBBtXmDsq;

		const SharedSpace *SharedBBtX = x->ObtainReadTempData("BBtX");
		const double *BBtX = SharedBBtX->ObtainReadData();
		double *temp = new double[r * r];

		// temp = etaxTV^T * BBtX, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transt, transn, &R, &R, &P, &one, const_cast<double *> (etaxTV), &P, const_cast<double *> (BBtX), &P, &zero, temp, &R);
		for (integer i = 0; i < r; i++)
		{
			temp[i + i * r] *= 2.0;
			for (integer j = i + 1; j < r; j++)
			{
				temp[i + j * r] += temp[j + i * r];
				temp[j + i * r] = temp[i + j * r];
			}
		}
		// exixTV = exixTV + 4 * mu * BBtX * temp, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &P, &R, &R, &fourmu, const_cast<double *> (BBtX), &P, temp, &R, &one, exixTV, &P);
		delete[] temp;
	};
}; /*end of ROPTLIB namespace*/
