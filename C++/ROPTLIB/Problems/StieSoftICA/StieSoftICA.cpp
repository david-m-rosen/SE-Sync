
#include "Problems/StieSoftICA/StieSoftICA.h"

/*Define the namespace*/
namespace ROPTLIB{

	StieSoftICA::StieSoftICA(double *inCs, integer inn, integer inp, integer inN)
	{
		Cs = inCs;
		n = inn;
		p = inp;
		N = inN;
	};

	double StieSoftICA::f(Variable *x) const
	{
		const double *xxM = x->ObtainReadData();
		SharedSpace *SharedCY = new SharedSpace(1, n * p * N);
		SharedSpace *SharedD = new SharedSpace(1, p * N);
		double *CY = SharedCY->ObtainWriteEntireData();
		double *D = SharedD->ObtainWriteEntireData();

		///*compute CY*/
		///*permute newCs = permute(Cs, [1, 3, 2])*/
		//double *newCs = new double[n * N * n + n * N * p];
		//double *newCY = newCs + n * N * n;
		//for (integer i = 0; i < n; i++)
		//{
		//	for (integer j = 0; j < N; j++)
		//	{
		//		dcopy_(&n, Cs + i * n + j * n * n, &GLOBAL::IONE, newCs + j * n + i * n * N, &GLOBAL::IONE);
		//	}
		//}

		//// newCs * Y
		//integer nN = n * N;
		//dgemm_(GLOBAL::N, GLOBAL::N, &nN, &p, &n, &GLOBAL::DONE, newCs, &nN, const_cast<double *> (xxM), &n, &GLOBAL::DZERO, newCY, &nN);

		///*permute CY = permute(newCY, [1, 3, 2])*/
		//for (integer i = 0; i < p; i++)
		//{
		//	for (integer j = 0; j < N; j++)
		//	{
		//		dcopy_(&n, newCY + j * n + i * N * n, &GLOBAL::IONE, CY + (i + j * p) * n, &GLOBAL::IONE);
		//	}
		//}
		//delete[] newCs;

		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		integer inc = 1;
		for (integer i = 0; i < N; i++)
		{
			// CY(0:n-1, 0:p-1, i) <- Cs(0:n-1, 0:n-1, i) * xxM, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
			dgemm_(transn, transn, &n, &p, &n, &one, Cs + n * n * i, &n, const_cast<double *> (xxM), &n, &zero, CY + n * p * i, &n);
		}

		for (integer i = 0; i < N; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				// output xxM(:, j)^T CY(:, j, i), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
				D[j + i * p] = ddot_(&n, const_cast<double *> (xxM + j * n), &GLOBAL::IONE, CY + n * p * i + j * n, &GLOBAL::IONE);
			}
		}
		integer pN = p * N;
		// output D(:)^T D(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		double result = -ddot_(&pN, D, &GLOBAL::IONE, D, &GLOBAL::IONE);
		if (UseGrad)
		{
			x->AddToTempData("CY", SharedCY);
			x->AddToTempData("D", SharedD);
		}
		else
		{
			delete SharedCY;
			delete SharedD;
		}
		return result;
	};

	void StieSoftICA::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *SharedCY = x->ObtainReadTempData("CY");
		const double *CY = SharedCY->ObtainReadData();
		const SharedSpace *SharedD = x->ObtainReadTempData("D");
		const double *D = SharedD->ObtainReadData();

		double *egfTV = egf->ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
		{
			egfTV[i] = 0;
		}
		integer np = n * p, inc = 1;
		double coef = 0;
		for (integer i = 0; i < N; i++)
		{
			for (integer j = 0; j < p; j++)
			{
				coef = -D[j + i * p] * 4;
				// egfTV(:, j) <- coef * CY(:, j, i) + egfTV(:, j), details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
				daxpy_(&n, &coef, const_cast<double *> (CY + i * n * p + j * n), &inc, egfTV + j * n, &inc);
			}
		}
	};

	//double StieSoftICA::f(Variable *x) const
	//{
	//	const double *xxM = x->ObtainReadData();
	//	SharedSpace *SharedCY = new SharedSpace(1, n * p * N);
	//	SharedSpace *SharedD = new SharedSpace(1, p * N);
	//	double *CY = SharedCY->ObtainWriteEntireData();
	//	double *D = SharedD->ObtainWriteEntireData();

	//	char *transn = const_cast<char *> ("n");
	//	double one = 1, zero = 0;
	//	integer inc = 1;
	//	for (integer i = 0; i < N; i++)
	//	{
	//		// CY(0:n-1, 0:p-1, i) <- Cs(0:n-1, 0:n-1, i) * xxM, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//		dgemm_(transn, transn, &n, &p, &n, &one, Cs + n * n * i, &n, const_cast<double *> (xxM), &n, &zero, CY + n * p * i, &n);
	//	}

	//	for (integer i = 0; i < N; i++)
	//	{
	//		for (integer j = 0; j < p; j++)
	//		{
	//			// output xxM(:, j)^T CY(:, j, i), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
	//			D[j + i * p] = ddot_(&n, const_cast<double *> (xxM + j * n), &inc, CY + n * p * i + j * n, &inc);
	//		}
	//	}
	//	integer pN = p * N;
	//	// output D(:)^T D(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
	//	double result = -ddot_(&pN, D, &inc, D, &inc);
	//	if (UseGrad)
	//	{
	//		x->AddToTempData("CY", SharedCY);
	//		x->AddToTempData("D", SharedD);
	//	}
	//	else
	//	{
	//		delete SharedCY;
	//		delete SharedD;
	//	}
	//	return result;
	//};

	//void StieSoftICA::EucGrad(Variable *x, Vector *egf) const
	//{
	//	const SharedSpace *SharedCY = x->ObtainReadTempData("CY");
	//	const double *CY = SharedCY->ObtainReadData();
	//	const SharedSpace *SharedD = x->ObtainReadTempData("D");
	//	const double *D = SharedD->ObtainReadData();

	//	double *egfTV = egf->ObtainWriteEntireData();
	//	for (integer i = 0; i < n * p; i++)
	//	{
	//		egfTV[i] = 0;
	//	}
	//	integer np = n * p, inc = 1;
	//	double coef = 0;
	//	for (integer i = 0; i < N; i++)
	//	{
	//		for (integer j = 0; j < p; j++)
	//		{
	//			coef = -D[j + i * p] * 4;
	//			// egfTV(:, j) <- coef * CY(:, j, i) + egfTV(:, j), details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
	//			daxpy_(&n, &coef, const_cast<double *> (CY + i * n * p + j * n), &inc, egfTV + j * n, &inc);
	//		}
	//	}
	//};

	void StieSoftICA::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const double *xxM = x->ObtainReadData();
		const double *etaxTV = etax->ObtainReadData();
		double *exixTV = exix->ObtainWriteEntireData();
		for (integer i = 0; i < n * p; i++)
		{
			exixTV[i] = 0;
		}

		const SharedSpace *SharedCY = x->ObtainReadTempData("CY");
		const double *CY = SharedCY->ObtainReadData();
		const SharedSpace *SharedD = x->ObtainReadTempData("D");
		const double *D = SharedD->ObtainReadData();

		double *temp = new double[n * p];
		double coef = 0;
		integer np = n * p, inc = 1;
		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		for (integer i = 0; i < N; i++)
		{
			// temp <- etaxTV, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
			dcopy_(&np, const_cast<double *> (etaxTV), &inc, temp, &inc);
			for (integer j = 0; j < p; j++)
			{
				// temp(:, j) <- D(j, i) * temp(:, j), details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
				dscal_(&n, const_cast<double *> (D + i * p + j), temp + j * n, &inc);
			}
			for (integer j = 0; j < p; j++)
			{
				// output etaxTV(:, j)^T CY(:, j, i), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
				coef = ddot_(&n, const_cast<double *> (etaxTV + j * n), &inc, const_cast<double *> (CY + i * n * p + j * n), &inc) * 2;
				// temp(:, j) <- coef * xxM(:, j) + temp(:, j), details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
				daxpy_(&n, &coef, const_cast<double *> (xxM + j * n), &inc, temp + j * n, &inc);
			}

			// exixTV <- exixTV + Cs(:, :, i) * temp, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
			dgemm_(transn, transn, &n, &p, &n, &one, Cs + i * n * n, &n, temp, &n, &one, exixTV, &n);
		}

		delete[] temp;

		Domain->ScaleTimesVector(x, -4.0, exix, exix);
	};
}; /*end of ROPTLIB namespace*/
