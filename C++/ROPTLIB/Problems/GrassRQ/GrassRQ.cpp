
#include "Problems/GrassRQ/GrassRQ.h"

/*Define the namespace*/
namespace ROPTLIB{

	GrassRQ::GrassRQ(double *inB, integer inn, integer inp)
	{
		B = inB;
		n = inn;
		p = inp;
	};

	GrassRQ::~GrassRQ(void)
	{
	};

	double GrassRQ::f(Variable *x) const
	{
		const double *xxM = x->ObtainReadData();
		Vector *Bx = x->ConstructEmpty();
		SharedSpace *Temp = new SharedSpace(Bx);
		double *temp = Bx->ObtainWriteEntireData();
		double result = 0;

		Matrix MB(B, n, n), MxxM(xxM, n, p), Mtemp(temp, n, p);
		// temp = B * xxM, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		Matrix::DGEMM(1, MB, false, MxxM, false, 0, Mtemp);

		integer length = n * p;
		// output temp(:)^T * xxM(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		result = ddot_(&length, temp, &GLOBAL::IONE, const_cast<double *> (xxM), &GLOBAL::IONE);
		if (UseGrad)
		{
			x->AddToTempData("Bx", Temp);
		}
		else
		{
			delete Temp;
		}
		return result;
	};

	void GrassRQ::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *Temp = x->ObtainReadTempData("Bx");
		Vector *Bx = Temp->GetSharedElement();
		Domain->ScaleTimesVector(x, 2.0, Bx, egf);
	};

	void GrassRQ::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const double *etaxTV = etax->ObtainReadData();
		double *exixTV = exix->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n");
		integer N = n, P = p, inc = 1, Length = N * P;
		double one = 1, zero = 0, two = 2;
		// exxiTV <- B * etaxTV, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &N, &P, &N, &one, B, &N, const_cast<double *> (etaxTV), &N, &zero, exixTV, &N);
		Domain->ScaleTimesVector(x, 2.0, exix, exix);
	};
}; /*end of ROPTLIB namespace*/
