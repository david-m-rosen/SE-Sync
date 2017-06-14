
#include "Problems/SPDTensorDL/SPDTensorDL.h"

/*Define the namespace*/
namespace ROPTLIB{

	SPDTensorDL::SPDTensorDL(double *inLs, double *inalpha, integer indim, integer inN, integer innum, double inlambdaX)
	{
		Ls = inLs;
		alpha = inalpha;
		dim = indim;
		N = inN;
		num = innum;
		lambdaX = inlambdaX;
	};

	SPDTensorDL::~SPDTensorDL(void)
	{
	};

	double SPDTensorDL::f(Variable *x) const
	{
		const double *xptr = x->ObtainReadData();
		/*each slice of Xaj is B alpha_j in [(6), CS15]*/
		SharedSpace *Xalpha = new SharedSpace(3, dim, dim, N);
		double *Xalphaptr = Xalpha->ObtainWriteEntireData();

		integer dd = dim * dim, nnum = num, NN = N;
		/*Xalpha <-- \mathbb{B} alpha*/
		dgemm_(GLOBAL::N, GLOBAL::N, &dd, &NN, &nnum, &GLOBAL::DONE, const_cast<double *> (xptr), &dd,
			alpha, &nnum, &GLOBAL::DZERO, Xalphaptr, &dd);
		x->AddToTempData("Xalpha", Xalpha);
		/*compute cholesky decomposition for all slices in Xalpha*/
		SPDTensor *Mani = dynamic_cast<SPDTensor *> (Domain);
		Mani->CholeskyRepresentation(x);

		const SharedSpace *SharedL = x->ObtainReadTempData("XaL");
		const double *L = SharedL->ObtainReadData();

		SharedSpace *SharedlogLXL = new SharedSpace(3, dim, dim, N);
		double *logLXL = SharedlogLXL->ObtainWriteEntireData();
		double *Ltmp = new double[dim * dim];
		integer length = dim * dim, ddim = dim, info;
		for (integer i = 0; i < N; i++)
		{
			dcopy_(&length, const_cast<double *> (L) + i * length, &GLOBAL::IONE, Ltmp, &GLOBAL::IONE);
			/*Solve the linear system Ls X = Li, i.e., X = Ls^{-1} Li. The solution X is stored in Li.
			Note that Li is a lower triangular matrix.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &ddim, &ddim, Ls + dim * dim * i, &ddim, Ltmp, &ddim, &info);
			if (info != 0)
			{
				printf("Warning: Solving linear system in SPDTensorDL::f failed with info:%d!\n", info);
			}
			dgemm_(GLOBAL::N, GLOBAL::T, &ddim, &ddim, &ddim, &GLOBAL::DONE, Ltmp, &ddim, Ltmp, &ddim, &GLOBAL::DZERO, logLXL + ddim * ddim * i, &ddim);
			Matrix MMt(logLXL + ddim * ddim * i, ddim, ddim);
			Matrix::LogSymmetricM(GLOBAL::L, MMt, MMt);
		}
		delete[] Ltmp;

		length = dim * dim * N;
		double result = dnrm2_(&length, logLXL, &GLOBAL::IONE);
		x->AddToTempData("logLXL", SharedlogLXL);
		result *= result;
		result /= 2.0;
		/*add \Omega(X) = \sum \tr(X_i)*/
		for (integer i = 0; i < num; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				result += lambdaX * xptr[i * dim * dim + j * dim + j];
			}
		}
		return result;
	};

	void SPDTensorDL::EucGrad(Variable *x, Vector *gf) const
	{
		const SharedSpace *SharedlogLXL = x->ObtainReadTempData("logLXL");
		const double *logLXL = SharedlogLXL->ObtainReadData();
		double *Log_Ainv_X_Xinv = new double[dim * dim * N];

		integer ddim = dim, info;

		const SharedSpace *SharedL = x->ObtainReadTempData("XaL");
		const double *Lx = SharedL->ObtainReadData();
		for (integer i = 0; i < N; i++)
		{
			/*tmp <-- log(Li^{-1} Xi Li^{-T}) Li^T */
			dgemm_(GLOBAL::N, GLOBAL::T, &ddim, &ddim, &ddim, &GLOBAL::DONE, const_cast<double *> (logLXL + ddim * ddim * i), &ddim,
				Ls + dim * dim * i, &ddim, &GLOBAL::DZERO, Log_Ainv_X_Xinv + dim * dim * i, &ddim);

			/*Solve the linear system Li^T X = tmp, i.e., X = Li^{-T} log(Li^{-1} X Li^{-T}) Li^T. The solution X is stored in tmp.
			Note that Li is a lower triangular matrix.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &ddim, &ddim, Ls + dim * dim * i, &ddim, Log_Ainv_X_Xinv + dim * dim * i, &ddim, &info);
			if (info != 0)
			{
				printf("The cholesky decompsotion in SPDTensorDL::EucGrad failed with info:%d!\n", info);
			}

			for (integer j = 0; j < dim; j++)
			{
				double swaptmp = 0;
				for (integer k = j; k < dim; k++)
				{
					swaptmp = (Log_Ainv_X_Xinv + dim * dim * i)[k + j * dim];
					(Log_Ainv_X_Xinv + dim * dim * i)[k + j * dim] = (Log_Ainv_X_Xinv + dim * dim * i)[j + k * dim];
					(Log_Ainv_X_Xinv + dim * dim * i)[j + k * dim] = swaptmp;
				}
			}

			/*Solve the linear system X Lx^T = tmp, i.e., X = Li^{-T} log(Li^{-1} X Li^{-T}) Li^T Lx^{-T}. The solution X is stored in tmp.
			We solve Lx X^t = tmp^T instead.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::N, GLOBAL::N, &ddim, &ddim, const_cast<double *> (Lx)+dim * dim * i, &ddim, Log_Ainv_X_Xinv + dim * dim * i, &ddim, &info);

			/*Solve the linear system X Lx = tmp, i.e., X = Li^{-T} log(Li^{-1} X Li^{-T}) Li^T Lx^{-T} Lx^{-1}. The solution X is stored in tmp.
			We can solve system Lx^T X^T = tmp^T. Since the Euclidean gradient is symmetric, we can solve Lx^T X = tmp^T instead.
			Details: http://www.netlib.org/lapack/explore-html/d6/d6f/dtrtrs_8f.html */
			dtrtrs_(GLOBAL::L, GLOBAL::T, GLOBAL::N, &ddim, &ddim, const_cast<double *> (Lx)+dim * dim * i, &ddim, Log_Ainv_X_Xinv + dim * dim * i, &ddim, &info);
		}
		double *gfVT = gf->ObtainWriteEntireData();
		
		integer dd = dim * dim, nnum = num, NN = N;
		dgemm_(GLOBAL::N, GLOBAL::T, &dd, &nnum, &NN, &GLOBAL::DONE, Log_Ainv_X_Xinv, &dd, alpha, &nnum, &GLOBAL::DZERO, gfVT, &dd);
		delete[] Log_Ainv_X_Xinv;

		/*add \nabla (\Omega(X))_i = I*/
		for (integer i = 0; i < num; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				gfVT[i * dim * dim + j * dim + j] += lambdaX;
			}
		}
	};

	void SPDTensorDL::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		printf("warning: SPDTensorDL::RieHessianEta has not been implemented!\n");
		etax->CopyTo(xix);
	};
}; /*end of ROPTLIB namespace*/
