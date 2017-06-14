
#include "Problems/EucPosSpCd/EucPosSpCd.h"

/*Define the namespace*/
namespace ROPTLIB{

	EucPosSpCd::EucPosSpCd(double *inLs, double *inB, double inlambdaa, integer indim, integer innum, integer inN)
	{
		Ls = inLs;
		B = inB;
		lambdaa = inlambdaa;
		dim = indim;
		num = innum;
		N = inN;
	};

	EucPosSpCd::~EucPosSpCd(void)
	{
	};

	double EucPosSpCd::f(Variable *x) const
	{
		const double *xptr = x->ObtainReadData();
		for (integer j = 0; j < N; j++)
		{
			double sumx = 0;
			for (integer i = 0; i < num; i++)
				sumx += xptr[i + num * j];
			if (sumx == 0)
				return 10000;
		}

		integer dd = dim * dim, nnum = num, ddim = dim, NN = N, info;
		SharedSpace *SharedBx = new SharedSpace(3, dim, dim, N);
		double *Bx = SharedBx->ObtainWriteEntireData();
		/*Bx <-- \sum_i B_i * x_i */
		dgemm_(GLOBAL::N, GLOBAL::N, &dd, &NN, &nnum, &GLOBAL::DONE, B, &dd, const_cast<double *> (xptr),
			&nnum, &GLOBAL::DZERO, Bx, &dd);

		/*Compute Cholesky decomposition for Bx*/
		SharedSpace *SharedBxL = new SharedSpace(3, dim, dim, N);
		double *BxL = SharedBxL->ObtainWriteEntireData();
		for (integer k = 0; k < N; k++)
		{
			/*Compute cholesky decomposition for k-th slice of Bx*/
			for (integer i = 0; i < dim; i++)
			{
				for (integer j = i; j < dim; j++)
				{
					BxL[i + j * dim + k * dim * dim] = 0;
					BxL[j + i * dim + k * dim * dim] = Bx[j + i * dim + k * dim * dim];
				}
			}
			dpotrf_(GLOBAL::L, &ddim, BxL + k * dim * dim, &ddim, &info);
			if (info != 0)
			{
				printf("Warning: Cholesky decompsition in EucPosSpCd::f failed with info:%d!\n", info);
			}
		}

		/*Compute Log(LXL)*/
		SharedSpace *SharedlogLXL = new SharedSpace(3, dim, dim, N);
		double *logLXL = SharedlogLXL->ObtainWriteEntireData();
		double *Ltmp = new double[dim * dim];
		for (integer k = 0; k < N; k++)
		{
			/*Solve the linear system L X = BxL, i.e., X = L^{-1} BxL. The solution X is stored in Ltmp.
			Note that Li is a lower triangular matrix.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dcopy_(&dd, BxL + k * dim * dim, &GLOBAL::IONE, Ltmp, &GLOBAL::IONE);
			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &ddim, &ddim, Ls + k * dim * dim, &ddim, Ltmp, &ddim, &info);
			if (info != 0)
			{
				printf("Warning: solving linear system in EucPosSpCd::f failed with info:%d!\n", info);
			}
			dgemm_(GLOBAL::N, GLOBAL::T, &ddim, &ddim, &ddim, &GLOBAL::DONE, Ltmp, &ddim, Ltmp, &ddim, &GLOBAL::DZERO, logLXL + k * dim * dim, &ddim);
			Matrix MMt(logLXL + k * dim * dim, ddim, ddim);
			Matrix::LogSymmetricM(GLOBAL::L, MMt, MMt);
		}

		delete[] Ltmp;
		
		integer length = dim * dim * N;
		double result = dnrm2_(&length, logLXL, &GLOBAL::IONE);
		result *= result;
		result /= 2.0;

		for (integer i = 0; i < num * N; i++)
			result += lambdaa * xptr[i];

		x->AddToTempData("logLXL", SharedlogLXL);
		x->AddToTempData("Bx", SharedBx); // this is not necessary
		x->AddToTempData("BxL", SharedBxL);
		return result;
	};

	void EucPosSpCd::EucGrad(Variable *x, Vector *gf) const
	{
		const SharedSpace *SharedlogLXL = x->ObtainReadTempData("logLXL");
		const double *logLXL = SharedlogLXL->ObtainReadData();
		const SharedSpace *SharedBxL = x->ObtainReadTempData("BxL");
		const double *BxL = SharedBxL->ObtainReadData();
		double *Log_Ainv_X_Xinv = new double[dim * dim * N];

		integer ddim = dim, info;

		for (integer i = 0; i < N; i++)
		{
			/*tmp <-- log(L^{-1} Bx L^{-T}) L^T */
			dgemm_(GLOBAL::N, GLOBAL::T, &ddim, &ddim, &ddim, &GLOBAL::DONE, const_cast<double *> (logLXL + i * dim * dim), &ddim,
				Ls + dim * dim * i, &ddim, &GLOBAL::DZERO, Log_Ainv_X_Xinv + i * dim * dim, &ddim);

			/*Solve the linear system L^T X = tmp, i.e., X = L^{-T} log(L^{-1} X L^{-T}) L^T. The solution X is stored in tmp.
			Note that L is a lower triangular matrix.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &ddim, &ddim, Ls + dim * dim * i, &ddim, Log_Ainv_X_Xinv + i * dim * dim, &ddim, &info);
			if (info != 0)
			{
				printf("The cholesky decompsotion in EucPosSpCd::EucGrad failed with info:%d!\n", info);
			}

			/*compute transpose of tmp*/
			for (integer j = 0; j < dim; j++)
			{
				double swaptmp = 0;
				for (integer k = j; k < dim; k++)
				{
					swaptmp = (Log_Ainv_X_Xinv)[k + j * dim + i * dim * dim];
					(Log_Ainv_X_Xinv)[k + j * dim + i * dim * dim] = (Log_Ainv_X_Xinv)[j + k * dim + i * dim * dim];
					(Log_Ainv_X_Xinv)[j + k * dim + i * dim * dim] = swaptmp;
				}
			}

			/*Solve the linear system X Lx^T = tmp, i.e., X = L^{-T} log(L^{-1} X L^{-T}) L^T Lx^{-T}. The solution X is stored in tmp.
			We solve Lx X^t = tmp^T instead.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &ddim, &ddim, const_cast<double *> (BxL + i * dim * dim), &ddim, Log_Ainv_X_Xinv + i * dim * dim, &ddim, &info);

			/*Solve the linear system X Lx = tmp, i.e., X = L^{-T} log(L^{-1} X L^{-T}) L^T Lx^{-T} Lx^{-1}. The solution X is stored in tmp.
			We can solve system Lx^T X^T = tmp^T. Since the Euclidean gradient is symmetric, we can solve Lx^T X = tmp^T instead.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &ddim, &ddim, const_cast<double *> (BxL + i * dim * dim), &ddim, Log_Ainv_X_Xinv + i * dim * dim, &ddim, &info);
		}

		double *gfTV = gf->ObtainWriteEntireData();
		integer dd = dim * dim;
		for (integer j = 0; j < N; j++)
		{
			for (integer i = 0; i < num; i++)
			{
				gfTV[i + j * num] = ddot_(&dd, Log_Ainv_X_Xinv + j * dim * dim, &GLOBAL::IONE, B + dim * dim * i, &GLOBAL::IONE) + lambdaa;
			}
		}

		delete[] Log_Ainv_X_Xinv;
	};

	void EucPosSpCd::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		printf("warning: EucPosSpCd::RieHessianEta has not been implemented!\n");
		etax->CopyTo(xix);
	};
}; /*end of ROPTLIB namespace*/
