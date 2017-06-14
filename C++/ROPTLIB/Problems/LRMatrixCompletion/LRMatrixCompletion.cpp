
#include "Problems/LRMatrixCompletion/LRMatrixCompletion.h"

/*Define the namespace*/
namespace ROPTLIB{

	LRMatrixCompletion::LRMatrixCompletion(integer *inir, integer *injc, double *inV, integer innz, integer inm, integer inn, integer inr)
	{
		ir = inir;
		jc = injc;
		V = inV;
		nz = innz;
		m = inm;
		n = inn;
		r = inr;
	};

	LRMatrixCompletion::~LRMatrixCompletion(void)
	{
	};

	void LRMatrixCompletion::ProjecOmegaUDVT(const double *U, const double *D, const double *V, integer inm, integer inn, integer inr, integer *inir, integer *injc, integer nz, double *result)
	{
		double *UD = new double[inm * inr];
		dgemm_(GLOBAL::N, GLOBAL::N, &inm, &inr, &inr, &GLOBAL::DONE, const_cast<double *> (U), &inm, const_cast<double *> (D), &inr, &GLOBAL::DZERO, UD, &inm);
		for (integer i = 0; i < nz; i++)
		{
			result[i] = 0;
			for (integer j = 0; j < inr; j++)
			{
				/*row: inir[i], col: injc[i]*/
				result[i] += UD[static_cast<integer> (inir[i]) + j * inm] * V[static_cast<integer> (injc[i]) + j * inn];
			}
		}
		delete[] UD;
	};

	double LRMatrixCompletion::f(Variable *x) const
	{
		ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
		const double *Uptr = ProdxxM->GetElement(0)->ObtainReadData();
		const double *Dptr = ProdxxM->GetElement(1)->ObtainReadData();
		const double *Vptr = ProdxxM->GetElement(2)->ObtainReadData();

		SharedSpace *EucRep = new SharedSpace(1, 2 + 3 * nz);
		double *EucRepptr = EucRep->ObtainWriteEntireData();
		EucRepptr[0] = 1; /*use sparse format*/
		EucRepptr[1] = static_cast<double> (nz); /*Set the number of nonzero entries*/
		ProjecOmegaUDVT(Uptr, Dptr, Vptr, m, n, r, ir, jc, nz, EucRepptr + 2);
		for (integer i = 0; i < nz; i++)
		{
			EucRepptr[2 + nz + i] = static_cast<double> (ir[i]);
			EucRepptr[2 + 2 * nz + i] = static_cast<double> (jc[i]);
		}
		daxpy_(const_cast<integer *> (&nz), &GLOBAL::DNONE, V, &GLOBAL::IONE, EucRepptr + 2, &GLOBAL::IONE);
		double result = 0.5 * ddot_(const_cast<integer *> (&nz), EucRepptr + 2, &GLOBAL::IONE, EucRepptr + 2, &GLOBAL::IONE);
		x->AddToTempData("EucRepinx", EucRep);
		return result;
	};
	
	void LRMatrixCompletion::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *EucRepinx = x->ObtainReadTempData("EucRepinx");
		SharedSpace *EucRep = EucRepinx->ConstructEmpty();
		EucRepinx->CopyTo(EucRep);
		egf->AddToTempData("EucRep", EucRep);
	};
	
	void LRMatrixCompletion::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		SharedSpace *EucRep = new SharedSpace(1, 2 + 3 * nz);
		double *EucRepptr = EucRep->ObtainWriteEntireData();
		EucRepptr[0] = 1; /*use sparse format*/
		EucRepptr[1] = static_cast<double> (nz); /*Set the number of nonzero entries*/
		for (integer i = 0; i < nz; i++)
		{
			EucRepptr[2 + nz + i] = static_cast<double> (ir[i]);
			EucRepptr[2 + 2 * nz + i] = static_cast<double> (jc[i]);
		}

		ProductElement *ProdxxM = dynamic_cast<ProductElement *>(x);
		const double *Uptr = ProdxxM->GetElement(0)->ObtainReadData();
		const double *Dptr = ProdxxM->GetElement(1)->ObtainReadData();
		const double *Vptr = ProdxxM->GetElement(2)->ObtainReadData();

		ProductElement *Prodetax = dynamic_cast<ProductElement *>(etax);
		const double *DUptr = Prodetax->GetElement(0)->ObtainReadData();
		const double *DDptr = Prodetax->GetElement(1)->ObtainReadData();
		const double *DVptr = Prodetax->GetElement(2)->ObtainReadData();

		ProjecOmegaUDVT(DUptr, Dptr, Vptr, m, n, r, ir, jc, nz, EucRepptr + 2);
		double *tmp = new double[nz];
		ProjecOmegaUDVT(Uptr, DDptr, Vptr, m, n, r, ir, jc, nz, tmp);
		daxpy_(const_cast<integer *> (&nz), &GLOBAL::DONE, tmp, &GLOBAL::IONE, EucRepptr + 2, &GLOBAL::IONE);
		ProjecOmegaUDVT(Uptr, Dptr, DVptr, m, n, r, ir, jc, nz, tmp);
		daxpy_(const_cast<integer *> (&nz), &GLOBAL::DONE, tmp, &GLOBAL::IONE, EucRepptr + 2, &GLOBAL::IONE);
		delete[] tmp;
		exix->AddToTempData("EucRep", EucRep);
	};
}; /*end of ROPTLIB namespace*/
