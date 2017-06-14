
#include "Problems/ObliqueSparsePCA/ObliqueSparsePCA.h"

/*Define the namespace*/
namespace ROPTLIB{

	ObliqueSparsePCA::ObliqueSparsePCA(double *inB, double *inDsq, double inmu, integer inp, integer inn, integer inr)
	{
		B = inB;
		Dsq = inDsq;
		p = inp;
		n = inn;
		r = inr;
		mu = inmu;
	};

	ObliqueSparsePCA::~ObliqueSparsePCA(void)
	{
	};

	double ObliqueSparsePCA::f(Variable *x) const
	{
		const double *xptr = x->ObtainReadData();

		double result = 0;
		//for (integer i = 0; i < p * r; i++)
		//{
		//	result += fabs(xptr[i]);
		//}
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

	void ObliqueSparsePCA::EucGrad(Variable *x, Vector *egf) const
	{
		const double *xptr = x->ObtainReadData();
		double *egfptr = egf->ObtainWriteEntireData();

		for (integer i = 0; i < p * r; i++)
			egfptr[i] = 0;

		//for (integer i = 0; i < p * r; i++)
		//{
		//	egfptr[i] = (xptr[i] > 0 ? 1 : -1);
		//	//if (xptr[i] > 1e-8)
		//	//	egfptr[i] = 1;
		//	//else
		//	//if (xptr[i] < -1e-8)
		//	//	egfptr[i] = -1;
		//	//else
		//	//{
		//	//	egfptr[i] = (genrandreal() < 0.5 ? 1 : -1);
		//	//}
		//}

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

	void ObliqueSparsePCA::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		printf("Warning: The cost function is not smooth. The action of the Hessian has not been implemented!\n");
		etax->CopyTo(exix);
	};
}; /*end of ROPTLIB namespace*/
