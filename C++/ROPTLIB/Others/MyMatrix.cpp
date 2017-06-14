
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	namespace GLOBAL{
		integer IZERO = 0, IONE = 1, ITWO = 2;
		double DZERO = 0, DONE = 1, DTWO = 2, DNONE = -1;
		doublecomplex ZZERO = { 0, 0 }, ZONE = { 1, 0 }, ZTWO = { 2, 0 }, ZNONE = { -1, 0 };
		char *N = const_cast<char *> ("n"), *T = const_cast<char *> ("t");
		char *L = const_cast<char *> ("l"), *R = const_cast<char *> ("r");
		char *V = const_cast<char *> ("v"), *C = const_cast<char *> ("c");
		char *U = const_cast<char *> ("u"), *A = const_cast<char *> ("a");
		char *S = const_cast<char *> ("s");
	};

	void Matrix::SPBtimesX(const double *B, const unsigned long long *ir, const unsigned long long *jc, integer nzmax, const double *X, integer m, integer n, integer p, double *result)
	{
		for (integer i = 0; i < m * p; i++)
			result[i] = 0;

		for (integer i = 0; i < n; i++)
		{
			for (unsigned long long j = jc[i]; j < jc[i + 1]; j++)
			{
				/*row: ir[j], column: i, entry: B[j]*/
				for (integer k = 0; k < p; k++)
				{
					result[ir[j] + k * m] += B[j] * X[i + k * n];
				}
			}
		}
	};

	void Matrix::XtimesSPB(const double *X, const double *B, const unsigned long long *ir, const unsigned long long *jc, integer nzmax, integer p, integer m, integer n, double *result)
	{
		for (integer i = 0; i < p * n; i++)
			result[i] = 0;

		for (integer i = 0; i < n; i++)
		{
			for (unsigned long long j = jc[i]; j < jc[i + 1]; j++)
			{
				/*row: ir[j], column: i, entry: B[j]*/
				for (integer k = 0; k < p; k++)
				{
					result[k + i * p] += X[k + ir[j] * p] * B[j];
				}
			}
		}
	};

	void Matrix::DGEMM(double alpha, Matrix &M1, bool trans1, Matrix &M2, bool trans2, double beta, Matrix &result)
	{
		if (trans1 && trans2)
		{
			if (M1.row != M2.col)
				printf("GEMM: the sizes of two matrices do not match!\n");
			dgemm_(GLOBAL::T, GLOBAL::T, &M1.col, &M2.row, &M1.row, &alpha, M1.matrix, &M1.inc, M2.matrix, &M2.inc, &beta, result.matrix, &result.inc);
		}
		else if (!trans1 && trans2)
		{
			if (M1.col != M2.col)
				printf("GEMM: the sizes of two matrices do not match!\n");
			dgemm_(GLOBAL::N, GLOBAL::T, &M1.row, &M2.row, &M1.col, &alpha, M1.matrix, &M1.inc, M2.matrix, &M2.inc, &beta, result.matrix, &result.inc);
		}
		else if (trans1 && !trans2)
		{
			if (M1.row != M2.row)
				printf("GEMM: the sizes of two matrices do not match!\n");
			dgemm_(GLOBAL::T, GLOBAL::N, &M1.col, &M2.col, &M1.row, &alpha, M1.matrix, &M1.inc, M2.matrix, &M2.inc, &beta, result.matrix, &result.inc);
		}
		else if (!trans1 && !trans2)
		{
			if (M1.col != M2.row)
				printf("GEMM: the sizes of two matrices do not match!\n");
			dgemm_(GLOBAL::N, GLOBAL::N, &M1.row, &M2.col, &M1.col, &alpha, M1.matrix, &M1.inc, M2.matrix, &M2.inc, &beta, result.matrix, &result.inc);
		}
		else
		{
			printf("impossible error!\n");
		}
	};

	void Matrix::CGEMM(doublecomplex alpha, Matrix &M1, bool trans1, Matrix &M2, bool trans2, doublecomplex beta, Matrix &result)
	{
		if (trans1 && trans2)
		{
			if (M1.row != M2.col)
				printf("GEMM: the sizes of two matrices do not match!\n");
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::C, GLOBAL::C, &M1.col, &M2.row, &M1.row, &alpha, (doublecomplex *)M1.matrix, &M1.inc,
				(doublecomplex *)M2.matrix, &M2.inc, &beta, (doublecomplex *)result.matrix, &result.inc);
#else
			zgemm_(GLOBAL::C, GLOBAL::C, &M1.col, &M2.row, &M1.row, (double *) &alpha, M1.matrix, &M1.inc,
				M2.matrix, &M2.inc, (double *) &beta, result.matrix, &result.inc);
#endif
		}
		else if (!trans1 && trans2)
		{
			if (M1.col != M2.col)
				printf("GEMM: the sizes of two matrices do not match!\n");
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::N, GLOBAL::C, &M1.row, &M2.row, &M1.col, &alpha, (doublecomplex *)M1.matrix, &M1.inc,
				(doublecomplex *)M2.matrix, &M2.inc, &beta, (doublecomplex *)result.matrix, &result.inc);
#else
			zgemm_(GLOBAL::N, GLOBAL::C, &M1.row, &M2.row, &M1.col, (double *) &alpha, M1.matrix, &M1.inc,
				M2.matrix, &M2.inc, (double *) &beta, result.matrix, &result.inc);
#endif
		}
		else if (trans1 && !trans2)
		{
			if (M1.row != M2.row)
				printf("GEMM: the sizes of two matrices do not match!\n");
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::C, GLOBAL::N, &M1.col, &M2.col, &M1.row, &alpha, (doublecomplex *)M1.matrix, &M1.inc,
				(doublecomplex *)M2.matrix, &M2.inc, &beta, (doublecomplex *)result.matrix, &result.inc);
#else
			zgemm_(GLOBAL::C, GLOBAL::N, &M1.col, &M2.col, &M1.row, (double *) &alpha, M1.matrix, &M1.inc,
				M2.matrix, &M2.inc, (double *) &beta, result.matrix, &result.inc);
#endif
		}
		else if (!trans1 && !trans2)
		{
			if (M1.col != M2.row)
				printf("GEMM: the sizes of two matrices do not match!\n");
#ifndef MATLAB_MEX_FILE
			zgemm_(GLOBAL::N, GLOBAL::N, &M1.row, &M2.col, &M1.col, &alpha, (doublecomplex *)M1.matrix, &M1.inc,
				(doublecomplex *)M2.matrix, &M2.inc, &beta, (doublecomplex *)result.matrix, &result.inc);
#else
			zgemm_(GLOBAL::N, GLOBAL::N, &M1.row, &M2.col, &M1.col, (double *) &alpha, M1.matrix, &M1.inc,
				M2.matrix, &M2.inc, (double *) &beta, result.matrix, &result.inc);
#endif
		}
		else
		{
			printf("impossible error!\n");
		}
	};

	void Matrix::CSYL(Matrix &A, Matrix &B, Matrix &C)
	{ // details can be found in http://www.netlib.org/lapack/lawnspdf/lawn75.pdf
		if (A.row != A.col || B.row != B.col)
			printf("CSYL: sizes do not match!\n");
		double *Bmatrix = nullptr;

		integer sdim, n = A.row, m = B.row;
		doublecomplex *eigsA = new doublecomplex[n + m + n * n + m * m];
		doublecomplex *eigsB = eigsA + n;
		doublecomplex *VSA = eigsB + m;
		doublecomplex *VSB = VSA + n * n;
		doublecomplex lworkopt;
		integer lwork = -1;
		integer info;
		double *rwork = new double[n];
		// compute the size of space required in the zgees
#ifndef MATLAB_MEX_FILE
		zgees_(GLOBAL::V, GLOBAL::N, nullptr, &n, (doublecomplex *)A.matrix, &A.inc, &sdim, eigsA, VSA, &n,
			&lworkopt, &lwork, rwork, nullptr, &info);
#else
		zgees_(GLOBAL::V, GLOBAL::N, nullptr, &n, A.matrix, &A.inc, &sdim, (double *) eigsA, (double *) VSA, &n,
			(double *) &lworkopt, &lwork, rwork, nullptr, &info);
#endif
		lwork = static_cast<integer> (lworkopt.r);
		doublecomplex *work = new doublecomplex[lwork];
		// generalized schur decomposition for matrices A.
		// details: http://www.netlib.org/lapack/explore-html/d8/d7e/zgees_8f.html
#ifndef MATLAB_MEX_FILE
		zgees_(GLOBAL::V, GLOBAL::N, nullptr, &n, (doublecomplex *)A.matrix, &A.inc, &sdim, eigsA, VSA, &n,
			work, &lwork, rwork, nullptr, &info);
#else
		zgees_(GLOBAL::V, GLOBAL::N, nullptr, &n, A.matrix, &A.inc, &sdim, (double *) eigsA, (double *) VSA, &n,
			(double *) work, &lwork, rwork, nullptr, &info);
#endif
		delete[] work;
		delete[] rwork;

		//printf("A:" << Matrix(A.matrix, 6, 2) << std::endl;
		//printf("VSA:" << Matrix((double *)VSA, 4, 2) << std::endl;

		if (B.matrix != A.matrix)
		{
			lwork = -1;
			rwork = new double[m];
			// compute the size of space required in the zgees
#ifndef MATLAB_MEX_FILE
			zgees_(GLOBAL::V, GLOBAL::N, nullptr, &m, (doublecomplex *)B.matrix, &B.inc, &sdim, eigsB, VSB, &m,
				&lworkopt, &lwork, rwork, nullptr, &info);
#else
			zgees_(GLOBAL::V, GLOBAL::N, nullptr, &m, B.matrix, &B.inc, &sdim, (double *) eigsB, (double *) VSB, &m,
				(double *) &lworkopt, &lwork, rwork, nullptr, &info);
#endif
			lwork = static_cast<integer> (lworkopt.r);
			work = new doublecomplex[lwork];
			// generalized schur decomposition for matrices B.
			// details: http://www.netlib.org/lapack/explore-html/d8/d7e/zgees_8f.html
#ifndef MATLAB_MEX_FILE
			zgees_(GLOBAL::V, GLOBAL::N, nullptr, &m, (doublecomplex *)B.matrix, &B.inc, &sdim, eigsB, VSB, &m,
				work, &lwork, rwork, nullptr, &info);
#else
			zgees_(GLOBAL::V, GLOBAL::N, nullptr, &m, B.matrix, &B.inc, &sdim, (double *) eigsB, (double *) VSB, &m,
				(double *) work, &lwork, rwork, nullptr, &info);
#endif
			delete[] work;
			delete[] rwork;
		}
		else
		{
			Bmatrix = new double[n * n * 2];
			integer length = 2 * n;
			for (integer i = 0; i < n; i++)
			{
				dcopy_(&length, A.matrix + A.inc * 2 * i, &GLOBAL::IONE, Bmatrix + 2 * n * i, &GLOBAL::IONE);
			}
			length = 2 * n * n;
			dscal_(&length, &GLOBAL::DNONE, Bmatrix, &GLOBAL::IONE);
			VSB = VSA;
		}

		//printf("C:" << Matrix(C.matrix, 6, 2) << std::endl;
		doublecomplex alpha = { 1, 0 }, beta = { 0, 0 };
		double *tempnbym = new double[2 * n * m];
		Matrix MVSA((double *)VSA, n, n), MC(C.matrix, n, m, C.inc), Mtempnbym(tempnbym, n, m);
		Matrix MVSB((double*)VSB, m, m);
		// tempnbym <- VSA^H * C
		CGEMM(alpha, MVSA, true, MC, false, beta, Mtempnbym);
		//printf("UA' * C:" << Matrix(tempnbym, 4, 2) << std::endl;
		// C <- tempnbym * VSB
		CGEMM(alpha, Mtempnbym, false, MVSB, false, beta, MC);

		//printf("UA' * C * UB:" << Matrix(C.matrix, 6, 2) << std::endl;

		doublecomplex *DIdentity = new doublecomplex[n * n + m * m + n * m];
		doublecomplex *EIdentity = DIdentity + n * n;
		doublecomplex *FZeros = EIdentity + m * m;
		integer length = 0;// 2 * (n * n + m * m + n * m);
		double scalar, dif;
		for (integer i = 0; i < n * n + m * m + n * m; i++)
		{
			DIdentity[i].r = 0;
			DIdentity[i].i = 0;
		}
		//dscal_(&length, &GLOBAL::DZERO, (double *)DIdentity, &GLOBAL::IONE);
		for (integer i = 0; i < n; i++)
			DIdentity[i + i * n].r = 1;
		for (integer i = 0; i < m; i++)
			EIdentity[i + i * m].r = 1;
		lwork = -1;
		integer *iwork = new integer[n + m + 2];
		if (A.matrix != B.matrix)
		{
			length = B.row * 2;
			for (integer i = 0; i < m; i++)
				dscal_(&length, &GLOBAL::DNONE, B.matrix + i * B.inc * 2, &GLOBAL::IONE);
			// compute the size of space required in the ztgsyl
#ifndef MATLAB_MEX_FILE
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, (doublecomplex *)A.matrix, &A.inc, (doublecomplex *)B.matrix, &B.inc, (doublecomplex *)C.matrix, &C.inc,
				DIdentity, &n, EIdentity, &m, FZeros, &n, &scalar, &dif, &lworkopt, &lwork, iwork, &info);
#else
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, A.matrix, &A.inc, B.matrix, &B.inc, C.matrix, &C.inc,
				(double *) DIdentity, &n, (double *) EIdentity, &m, (double *) FZeros, &n, &scalar, &dif, (double *) &lworkopt, &lwork, iwork, &info);
#endif
			lwork = static_cast<integer> (lworkopt.r);
			work = new doublecomplex[lwork];
			// generalized Sylvester equation.
			// details: http://www.netlib.org/lapack/explore-html/d8/d7e/zgees_8f.html
#ifndef MATLAB_MEX_FILE
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, (doublecomplex *)A.matrix, &A.inc, (doublecomplex *)B.matrix, &B.inc, (doublecomplex *)C.matrix, &C.inc,
				DIdentity, &n, EIdentity, &m, FZeros, &n, &scalar, &dif, work, &lwork, iwork, &info);
#else
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, A.matrix, &A.inc, B.matrix, &B.inc, C.matrix, &C.inc,
				(double *)DIdentity, &n, (double *)EIdentity, &m, (double *)FZeros, &n, &scalar, &dif, (double *)work, &lwork, iwork, &info);
#endif
			delete[] work;
		}
		else
		{
			// compute the size of space required in the ztgsyl
#ifndef MATLAB_MEX_FILE
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, (doublecomplex *)A.matrix, &A.inc, (doublecomplex *)Bmatrix, &m, (doublecomplex *)C.matrix, &C.inc,
				DIdentity, &n, EIdentity, &m, FZeros, &n, &scalar, &dif, &lworkopt, &lwork, iwork, &info);
#else
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, A.matrix, &A.inc, Bmatrix, &m, C.matrix, &C.inc,
				(double *)DIdentity, &n, (double *)EIdentity, &m, (double *)FZeros, &n, &scalar, &dif, (double *)&lworkopt, &lwork, iwork, &info);
#endif
			lwork = static_cast<integer> (lworkopt.r);
			work = new doublecomplex[lwork];
			// generalized Sylvester equation.
			// details: http://www.netlib.org/lapack/explore-html/d8/d7e/zgees_8f.html
#ifndef MATLAB_MEX_FILE
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, (doublecomplex *)A.matrix, &A.inc, (doublecomplex *)Bmatrix, &m, (doublecomplex *)C.matrix, &C.inc,
				DIdentity, &n, EIdentity, &m, FZeros, &n, &scalar, &dif, work, &lwork, iwork, &info);
#else
			ztgsyl_(GLOBAL::N, &GLOBAL::IZERO, &n, &m, A.matrix, &A.inc, Bmatrix, &m, C.matrix, &C.inc,
				(double *)DIdentity, &n, (double *)EIdentity, &m, (double *)FZeros, &n, &scalar, &dif, (double *)work, &lwork, iwork, &info);
#endif
			delete[] work;
		}
		delete[] iwork;

		// tempnbym <- VSA^H * C
		CGEMM(alpha, MVSA, false, MC, false, beta, Mtempnbym);
		// C <- tempnbym * VSB
		CGEMM(alpha, Mtempnbym, false, MVSB, true, beta, MC);

		delete[] DIdentity;
		delete[] eigsA;
		delete[] tempnbym;
		if (B.matrix != nullptr)
		{
			delete[] Bmatrix;
		}
	};


	void Matrix::EigenSymmetricM(char *UorL, Matrix &S, Matrix &eigenvalues, Matrix &eigenvectors)
	{
		/* eigenvalue decomposition: dsyev approach. */
		integer N = S.row, inc = S.inc;
		for (integer i = 0; i < N; i++)
			dcopy_(&N, S.matrix + i * inc, &GLOBAL::IONE, eigenvectors.matrix + i * eigenvectors.inc, &GLOBAL::IONE);
		integer lwork = -1, info;
		double lworkopt;
		dsyev_(GLOBAL::V, UorL, &N, eigenvectors.matrix, &eigenvectors.inc, eigenvalues.matrix, &lworkopt, &lwork, &info);
		/*allocate the desired memory*/
		lwork = static_cast<integer> (lworkopt);
		double *work = new double[lwork];
		dsyev_(GLOBAL::V, UorL, &N, eigenvectors.matrix, &eigenvectors.inc, eigenvalues.matrix, work, &lwork, &info);
		delete[] work;

		/////* eigenvalue decomposition: dsyevr approach, not numerically stable. */
		//integer N = S.row, inc = S.inc, numeig;
		//integer *isuppz = new integer[2 * N];
		//integer lwork = -1, liwork = -1, liworkopt, info;
		//double lworkopt;
		//dsyevr_(GLOBAL::V, GLOBAL::A, UorL, &N, S.matrix, &inc, &GLOBAL::DZERO, &GLOBAL::DZERO, &GLOBAL::IZERO, &GLOBAL::IZERO,
		//	&GLOBAL::DZERO, &numeig, eigenvalues.matrix, eigenvectors.matrix, &eigenvectors.inc, isuppz, 
		//	&lworkopt, &lwork, &liworkopt, &liwork, &info);
		//lwork = static_cast<integer> (lworkopt);
		//liwork = liworkopt;
		//double *work = new double[lwork];
		//integer *iwork = new integer[liwork];
		//dsyevr_(GLOBAL::V, GLOBAL::A, UorL, &N, S.matrix, &inc, &GLOBAL::DZERO, &GLOBAL::DZERO, &GLOBAL::IZERO, &GLOBAL::IZERO,
		//	&GLOBAL::DZERO, &numeig, eigenvalues.matrix, eigenvectors.matrix, &eigenvectors.inc, isuppz,
		//	work, &lwork, iwork, &liwork, &info);
		//if (info != 0)
		//	printf("Matrix::EigenSymmetricM failed!\n");
		//delete[] iwork;
		//delete[] work;
		//delete[] isuppz;


		///* eigenvalue decomposition: dsyevd approach, not numerically stable. */
		//integer N = S.row, inc = S.inc;
		//integer info;
		//integer lwork = -1, liwork = -1, liworkopt;
		//double lworkopt;
		//for (integer i = 0; i < N; i++)
		//	dcopy_(&N, S.matrix + i * inc, &GLOBAL::IONE, eigenvectors.matrix + i * eigenvectors.inc, &GLOBAL::IONE);
		//dsyevd_(GLOBAL::V, UorL, &N, eigenvectors.matrix, &eigenvectors.inc, eigenvalues.matrix, &lworkopt, &lwork, &liworkopt, &liwork, &info);
		//lwork = static_cast<integer> (lworkopt);
		//liwork = liworkopt;
		//double *work = new double[lwork];
		//integer *iwork = new integer[liwork];
		//dsyevd_(GLOBAL::V, UorL, &N, eigenvectors.matrix, &eigenvectors.inc, eigenvalues.matrix, work, &lwork, iwork, &liwork, &info);
		//if (info != 0)
		//	printf("Matrix::EigenSymmetricM failed!\n");
		//delete[] work;
		//delete[] iwork;

		///*Schur decomposition approach:*/
		//integer N = S.row, inc = S.inc, sdim;
		//double *wi = new double[N];
		//integer info;
		//integer lwork = -1;
		//double lworkopt;
		///*Using lwork = -1 to find the optimal required space.*/
		//dgees_(GLOBAL::V, GLOBAL::N, nullptr, &N, S.matrix, &inc, &sdim, eigenvalues.matrix, wi, eigenvectors.matrix, &eigenvectors.inc,
		//	&lworkopt, &lwork, nullptr, &info);
		///*allocate the desired memory*/
		//lwork = static_cast<integer> (lworkopt);
		//double *work = new double[lwork];
		///*Compute the schur decomposition, S = eigenvectors * diag(eigenvalues) * eigenvectors^T.
		//Details: http://www.netlib.org/lapack/explore-html/d9/d8e/group__double_g_eeigen.html */
		//dgees_(GLOBAL::V, GLOBAL::N, nullptr, &N, S.matrix, &inc, &sdim, eigenvalues.matrix, wi, eigenvectors.matrix, &eigenvectors.inc,
		//	work, &lwork, nullptr, &info);
		//if (info != 0)
		//	printf("Matrix::EigenSymmetricM failed!\n");
		//delete[] work;
		//delete[] wi;

		///* eigenvalue decomposition: dsyevx_ approach */
		//integer N = S.row, inc = S.inc, numeigvalues;
		//double tol = dlamch_(GLOBAL::S);
		//integer info;
		//integer lwork = -1;
		//double lworkopt;
		//integer *iwork = new integer[6 * N];
		//integer *ifail = iwork + 5 * N;
		///*Using lwork = -1 to find the optimal required space.*/
		//dsyevx_(GLOBAL::V, GLOBAL::A, UorL, &N, S.matrix, &inc, &GLOBAL::DZERO, &GLOBAL::DZERO, &GLOBAL::IZERO, &GLOBAL::IZERO,
		//	&tol, &numeigvalues, eigenvalues.matrix, eigenvectors.matrix, &eigenvectors.inc, &lworkopt, &lwork, iwork, ifail, &info);
		///*allocate the desired memory*/
		//lwork = static_cast<integer> (lworkopt);
		//double *work = new double[lwork];
		///*Compute the eigenvalue decomposition, S = eigenvectors * diag(eigenvalues) * eigenvectors^T.
		//Details: http://www.netlib.org/lapack/explore-html/d2/d8a/group__double_s_yeigen.html*/
		//dsyevx_(GLOBAL::V, GLOBAL::A, UorL, &N, S.matrix, &inc, &GLOBAL::DZERO, &GLOBAL::DZERO, &GLOBAL::IZERO, &GLOBAL::IZERO,
		//	&tol, &numeigvalues, eigenvalues.matrix, eigenvectors.matrix, &eigenvectors.inc, work, &lwork, iwork, ifail, &info);

		//if (info != 0)
		//{
		//	printf("Warning: The lapack function dsyevx_ used in Matrix::EigenSymmetricM did not secceessfully converge!\n");
		//}

		//delete[] iwork;
		//delete[] work;
	};

	void Matrix::ExpSymmetricM(char *UorL, Matrix &S, Matrix &ExpS)
	{
		integer N = S.row;
		double *eigenvalues = new double[N + 2 * N * N];
		double *eigenvectors = eigenvalues + N;
		double *eigenvectorsD = eigenvectors + N * N;
		Matrix E(eigenvalues, N, 1), V(eigenvectors, N, N), VD(eigenvectorsD, N, N);
		EigenSymmetricM(UorL, S, E, V);
		integer length = N * N;
		dcopy_(&length, eigenvectors, &GLOBAL::IONE, VD.matrix, &GLOBAL::IONE);
		for (integer i = 0; i < N; i++)
		{
			double a = exp(eigenvalues[i]);
			dscal_(&N, &a, eigenvectors + i * N, &GLOBAL::IONE);
		}
		DGEMM(GLOBAL::DONE, V, false, VD, true, GLOBAL::DZERO, ExpS);

		delete[] eigenvalues;
	};

	void Matrix::LogSymmetricM(char *UorL, Matrix &S, Matrix &LogS)
	{
		integer N = S.row;
		double *eigenvalues = new double[N + 2 * N * N];
		double *eigenvectors = eigenvalues + N;
		double *eigenvectorsD = eigenvectors + N * N;
		Matrix E(eigenvalues, N, 1), V(eigenvectors, N, N), VD(eigenvectorsD, N, N);
		EigenSymmetricM(UorL, S, E, V);
		integer length = N * N;
		dcopy_(&length, eigenvectors, &GLOBAL::IONE, VD.matrix, &GLOBAL::IONE);
		for (integer i = 0; i < N; i++)
		{
			if (eigenvalues <= 0)
			{
				printf("Error: The matrix for Matrix::LogSymmetricM is not symmetric positive definite!!\n");
				return;
			}
			double a = log(eigenvalues[i]);
			dscal_(&N, &a, eigenvectors + i * N, &GLOBAL::IONE);
		}
		DGEMM(GLOBAL::DONE, V, false, VD, true, GLOBAL::DZERO, LogS);

		delete[] eigenvalues;
	};

	const Matrix &Matrix::operator=(const Matrix &right)
	{
		assert(row == right.row && col == right.col);

		for (int i = 0; i < row; i++)
			for (int j = 0; j < col; j++)
				matrix[i + j * inc] = right.matrix[i + j * right.inc];

		return *this;
	}
}; /*end of ROPTLIB namespace*/
