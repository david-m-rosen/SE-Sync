
#include "Problems/WeightedLowrank/WeightedLowRank.h"

/*Define the namespace*/
namespace ROPTLIB{

	WeightedLowRank::WeightedLowRank(double *inA, double *inW, integer inm, integer inn, integer inr)
	{
		A = inA;
		W = inW;
		m = inm;
		n = inn;
		r = inr;
	};

	WeightedLowRank::~WeightedLowRank(void)
	{
	};

	double WeightedLowRank::f(Variable *x) const
	{
		ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
		const double *Uptr = ProdxxM->GetElement(0)->ObtainReadData();
		const double *Dptr = ProdxxM->GetElement(1)->ObtainReadData();
		const double *Vptr = ProdxxM->GetElement(2)->ObtainReadData();

		char *transn = const_cast<char *> ("n");
		char *transt = const_cast<char *> ("t");
		char *uplo = const_cast<char *> ("u");
		double one = 1, zero = 0, neg_one = -1;
		integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r;

		double *UDptr = new double[MR];
		// UDptr <- Uptr * Dptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &M, &R, &R, &one, const_cast<double *> (Uptr), &M, const_cast<double *> (Dptr), &R, &zero, UDptr, &M);
		SharedSpace *Temp1 = new SharedSpace(2, m, n);
		double *Xptr = Temp1->ObtainWriteEntireData();
		// Xptr <- UDptr * Vptr^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transt, &M, &N, &R, &one, UDptr, &M, const_cast<double *> (Vptr), &N, &zero, Xptr, &M);
		delete[] UDptr;

		SharedSpace *Temp2 = new SharedSpace(2, m, n);
		double *Errptr = Temp2->ObtainWriteEntireData();
		// Errptr <- A, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&MN, A, &inc, Errptr, &inc);
		// Errptr <- Errptr - Xptr, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&MN, &neg_one, Xptr, &inc, Errptr, &inc);

		SharedSpace *Temp3 = new SharedSpace(2, m, n);
		double *QXptr = Temp3->ObtainWriteEntireData();
		// QXptr = W * Errptr, details: http://www.netlib.org/lapack/explore-html/d8/dbe/dsymv_8f.html
		dsymv_(uplo, &MN, &one, W, &MN, Errptr, &inc, &zero, QXptr, &inc);

		double result = 0;
		// compute Errptr(:)^T QXptr(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		result = ddot_(&MN, Errptr, &inc, QXptr, &inc);
		if (UseGrad)
		{
			x->AddToTempData("X", Temp1);
			x->AddToTempData("err", Temp2);
			x->AddToTempData("QX", Temp3);
		}
		else
		{
			delete Temp1;
			delete Temp2;
			delete Temp3;
		}
		return result;
	};
	
	void WeightedLowRank::EucGrad(Variable *x, Vector *egf) const
	{
		ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
	    double *Uptr = const_cast<double *> (ProdxxM->GetElement(0)->ObtainReadData());
	   	double *Dptr = const_cast<double *> (ProdxxM->GetElement(1)->ObtainReadData());
	    double *Vptr = const_cast<double *> (ProdxxM->GetElement(2)->ObtainReadData());
		char *transn = const_cast<char *> ("n");
	  	char *transt = const_cast<char *> ("t");
	    double one = 1, zero = 0, neg_one = -1.0, neg_two = -2.0;
	    integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r, NR = n * r; 
	
	    const SharedSpace *Temp = x->ObtainReadTempData("QX");
		const double *QXptr = Temp->ObtainReadData();

		SharedSpace *EucRep = new SharedSpace(1, 1 + m * n);
		double *EucRepptr = EucRep->ObtainWriteEntireData();
		EucRepptr[0] = 0;
		double *fullgrad = EucRepptr + 1;

		dcopy_(&MN, const_cast<double *> (QXptr), &inc, fullgrad, &inc);
		dscal_(&MN, &neg_two, fullgrad, &inc);

		egf->AddToTempData("EucRep", EucRep);
	};
	
	void WeightedLowRank::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		dynamic_cast<LowRank *> (Domain)->ExtrToEucRep(x, etax);
		const SharedSpace *EucRep = etax->ObtainReadTempData("EucRep");
		const double *EucRepptr = EucRep->ObtainReadData();

		SharedSpace *EucRep_R = new SharedSpace(1, 1 + m * n);
		double *EucRep_Rptr = EucRep_R->ObtainWriteEntireData();
		EucRep_Rptr[0] = 0;
		integer MN = n * m;
		// exix = 2 * W * etax, details: http://www.netlib.org/lapack/explore-html/d8/dbe/dsymv_8f.html
		dsymv_(GLOBAL::U, &MN, &GLOBAL::DTWO, W, &MN, const_cast<double *> (EucRepptr + 1), &GLOBAL::IONE, &GLOBAL::DZERO, EucRep_Rptr + 1, &GLOBAL::IONE);

		exix->AddToTempData("EucRep", EucRep_R);
	};

	//void WeightedLowRank::RieGrad(Variable *x, Vector *gf) const
	//{
	//	ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
	//	const double *Uptr = ProdxxM->GetElement(0)->ObtainReadData();
	//	const double *Dptr = ProdxxM->GetElement(1)->ObtainReadData();
	//	const double *Vptr = ProdxxM->GetElement(2)->ObtainReadData();
	//	char *transn = const_cast<char *> ("n");
	//	char *transt = const_cast<char *> ("t");
	//	double one = 1, zero = 0, neg_one = -1.0, neg_two = -2.0;
	//	integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r, NR = n * r, RR = r * r;
	//	const SharedSpace *Temp = x->ObtainReadTempData("QX");
	//	const double *QXptr = Temp->ObtainReadData();
	//	double *fullgrad = new double[MN];
	//	// fullgrad <- Qxptr, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
	//	dcopy_(&MN, const_cast<double *> (QXptr), &inc, fullgrad, &inc);
	//	// fullgrad <- -2 * fullgrad, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
	//	dscal_(&MN, &neg_two, fullgrad, &inc);
	//	double *XiVptr = new double[MR];
	//	// XiVptr <- fullgrad * Vptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//	dgemm_(transn, transn, &M, &R, &N, &one, fullgrad, &M, const_cast<double *> (Vptr), &N, &zero, XiVptr, &M);

	//	double *XiUptr = new double[NR];
	//	// XiUptr <- fullgrad^T * Uptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//	dgemm_(transt, transn, &N, &R, &M, &one, fullgrad, &M, const_cast<double *> (Uptr), &M, &zero, XiUptr, &N);
	//	delete[] fullgrad;

	//	// compute inverse of D
	//	integer info;
	//	integer *IPIV = new integer[R + 1];
	//	double *WORK = new double[RR];
	//	double *Dinv = new double[RR];
	//	// Dinv <- Dptr, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
	//	dcopy_(&RR, const_cast<double *> (Dptr), &inc, Dinv, &inc);
	//	// LU factorization for Dinv, LU factors are stored in Dinv,
	//	// details: http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
	//	dgetrf_(&R, &R, Dinv, &R, IPIV, &info);
	//	// Compute the inverse of Dinv based on Dinv's LU factorization,
	//	// details: http://www.netlib.org/lapack/explore-html/df/da4/dgetri_8f.html
	//	dgetri_(&R, Dinv, &R, IPIV, WORK, &RR, &info);
	//	delete[] IPIV;
	//	delete[] WORK;

	//	double *gfTV = gf->ObtainWriteEntireData();
	//	double *Udotptr = gfTV;
	//	double *Ddotptr = gfTV + m * r;
	//	double *Vdotptr = Ddotptr + r * r;
	//	// Ddotptr <- Uptr^T * XiVptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//	dgemm_(transt, transn, &R, &R, &M, &one, const_cast<double *> (Uptr), &M, XiVptr, &M, &zero, Ddotptr, &R);
	//	// Udotptr <- Uptr * Ddotptr, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//	dgemm_(transn, transn, &M, &R, &R, &one, const_cast<double *> (Uptr), &M, Ddotptr, &R, &zero, Udotptr, &M);
	//	// Udotptr <- -1 * Udotptr, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
	//	dscal_(&MR, &neg_one, Udotptr, &inc);
	//	// Udotptr <- XiVptr + Udotptr, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
	//	daxpy_(&MR, &one, XiVptr, &inc, Udotptr, &inc);
	//	// Vdotptr <- Vptr * Ddotptr^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//	dgemm_(transn, transt, &N, &R, &R, &one, const_cast<double *> (Vptr), &N, Ddotptr, &R, &zero, Vdotptr, &N);
	//	// Vdotptr <- -1 * Vdotptr, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
	//	dscal_(&NR, &neg_one, Vdotptr, &inc);
	//	// Vdotptr <- XiUptr + Vdotptr, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
	//	daxpy_(&NR, &one, XiUptr, &inc, Vdotptr, &inc);

	//	double *Udottemp = new double[MR];
	//	double *Vdottemp = new double[NR];
	//	// Udottemp <- Udotptr * Dinv, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//	dgemm_(transn, transn, &M, &R, &R, &one, Udotptr, &M, Dinv, &R, &zero, Udottemp, &M);
	//	// Vdottemp <- Vdotptr * Dinv^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
	//	dgemm_(transn, transt, &N, &R, &R, &one, Vdotptr, &N, Dinv, &R, &zero, Vdottemp, &N);
	//	// Udotptr <- Udottemp, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
	//	dcopy_(&MR, Udottemp, &inc, Udotptr, &inc);
	//	// Vdotptr <- Vdottemp, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
	//	dcopy_(&NR, Vdottemp, &inc, Vdotptr, &inc);
	//	//gf->Print("gf");//---
	//	delete[] Udottemp;
	//	delete[] Vdottemp;
	//	delete[] Dinv;
	//	delete[] XiUptr;
	//	delete[] XiVptr;
	//};

	//void WeightedLowRank::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	//{
	//	etax->CopyTo(xix);

	//	/*ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
	//	 double *Uptr = const_cast<double *> (ProdxxM->GetElement(0)->ObtainReadData());
	//	 double *Dptr = const_cast<double *> (ProdxxM->GetElement(1)->ObtainReadData());
	//	 double *Vptr = const_cast<double *> (ProdxxM->GetElement(2)->ObtainReadData());
	//	 char *transn = const_cast<char *> ("n");
	//	 char *transt = const_cast<char *> ("t");
	//	 double one = 1, zero = 0, neg_one = -1.0, neg_two = -2.0;
	//	 integer inc = 1, M = m, R = r, N = n, MN = m * n, MR = m * r, NR = n * r;

	//	 const SharedSpace *Temp = x->ObtainReadTempData("QX");
	//	 const double *QXptr = Temp->ObtainReadData();

	//	 double *gfTV = etax->ObtainWriteEntireData();
	//	 double *Udotptr = gfTV;
	//	 double *Ddotptr = gfTV + m * r;
	//	 double *Vdotptr = Ddotptr + r * r;

	//	 double *UDdot = new double[MR];
	//	 dgemm_(transn, transn, &M, &R, &R, &one, Uptr, &M, Ddotptr, &R, &zero, UDdot, &M);

	//	 double *Qeta = new double[MN];
	//	 dgemm_(transn, transt, &M, &N, &R, &one, Udotptr, &M, Vptr, &N, &zero, Qeta, &M);*/






	//};

}; /*end of ROPTLIB namespace*/
