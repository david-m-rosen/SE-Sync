
#include "Problems/StieSparseBrockett/StieSparseBrockett.h"

/*Define the namespace*/
namespace ROPTLIB{

	StieSparseBrockett::StieSparseBrockett(double *inB, unsigned long long *inir, unsigned long long *injc, integer innzmax, double *inD, integer inn, integer inp)
	{
		B = inB;
		ir = inir;
		jc = injc;
		nzmax = innzmax;
		D = inD;
		n = inn;
		p = inp;
	};

	StieSparseBrockett::~StieSparseBrockett(void)
	{
	};

	double StieSparseBrockett::f(Variable *x) const
	{
		const double *xxM = x->ObtainReadData();
		Vector *BxD = x->ConstructEmpty();
		SharedSpace *Temp = new SharedSpace(BxD);
		double *temp = BxD->ObtainWriteEntireData();
		double result = 0;

		/*Temp <-- B * X */
		Matrix::SPBtimesX(B, ir, jc, nzmax, xxM, n, n, p, temp);

		for (integer i = 0; i < p; i++)
		{
			// temp(i * n : i * n + N - 1) <- D[i] * temp(i * n : i * n + N - 1),
			// details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(const_cast<integer *> (&n), &D[i], temp + i * n, &GLOBAL::IONE);
		}
		integer length = n * p;
		// output temp(:)^T * xxM(:), details: http://www.netlib.org/lapack/explore-html/d5/df6/ddot_8f.html
		result = ddot_(&length, temp, &GLOBAL::IONE, const_cast<double *> (xxM), &GLOBAL::IONE);
		if (UseGrad)
		{
			x->AddToTempData("BxD", Temp);
		}
		else
		{
			delete Temp;
		}
		return result;
	};

	void StieSparseBrockett::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *Temp = x->ObtainReadTempData("BxD");
		Vector *BxD = Temp->GetSharedElement();
		Domain->ScaleTimesVector(x, 2.0, BxD, egf);
	};

	void StieSparseBrockett::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const double *etaxTV = etax->ObtainReadData();
		double *exixTV = exix->ObtainWriteEntireData();

		char *transn = const_cast<char *> ("n");
		integer N = n, P = p, inc = 1, Length = N * P;
		double one = 1, zero = 0, two = 2;

		/*exxiTV <- B * etaxTV */
		Matrix::SPBtimesX(B, ir, jc, nzmax, etaxTV, n, n, p, exixTV);

		for (integer i = 0; i < p; i++)
		{
			// exixTV(i * n : i * n + N - 1) <- D[i] * exixTV(i * n : i * n + N - 1),
			// details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
			dscal_(&N, &D[i], exixTV + i * n, &inc);
		}
		Domain->ScaleTimesVector(x, 2.0, exix, exix);
	};
}; /*end of ROPTLIB namespace*/
