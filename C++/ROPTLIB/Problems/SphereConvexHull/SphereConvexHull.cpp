
#include "Problems/SphereConvexHull/SphereConvexHull.h"

/*Define the namespace*/
namespace ROPTLIB{

	SphereConvexHull::SphereConvexHull(const Manifold *inMani, Variable *inx, Vector **inW, integer inlengthW, QuasiNewton *insolver, void (QuasiNewton::*inHv)(Vector *v, Vector *result))
	{
		Mani = inMani;
		x = inx;
		W = inW;
		lengthW = inlengthW;
		solver = insolver;
		Hv = inHv;
	};

	SphereConvexHull::~SphereConvexHull(void)
	{
	};

	double SphereConvexHull::f(Variable *x) const
	{
		Vector *xsq = x->ConstructEmpty();
		x->CopyTo(xsq);
		double *xsqPtr = xsq->ObtainWritePartialData();
		for (integer i = 0; i < xsq->Getlength(); i++)
			xsqPtr[i] *= xsqPtr[i];

		Vector *Wxsq = W[0]->ConstructEmpty();
		SharedSpace *SharedWxsq = new SharedSpace(Wxsq);
		Mani->ScaleTimesVector(x, xsqPtr[0], W[0], Wxsq);
		for (integer i = 1; i < lengthW; i++)
		{
			Mani->scalarVectorAddVector(x, xsqPtr[i], W[i], Wxsq, Wxsq);
		}
		delete xsq;
		Vector *PWxsq = Wxsq->ConstructEmpty();
		SharedSpace *SharedPWxsq = new SharedSpace(PWxsq);
		if (Hv == NULL)//-- nullptr)
			Wxsq->CopyTo(PWxsq);
		else
			(solver->*Hv)(Wxsq, PWxsq);

		const double *WxsqPtr = Wxsq->ObtainReadData();
		const double *PWxsqPtr = PWxsq->ObtainReadData();
		integer length = Wxsq->Getlength();
		double result = ddot_(&length, const_cast<double *> (WxsqPtr), &GLOBAL::IONE, const_cast<double *> (PWxsqPtr), &GLOBAL::IONE);

		x->AddToTempData("Wxsq", SharedWxsq);
		x->AddToTempData("PWxsq", SharedPWxsq);

		return result;
	};

	void SphereConvexHull::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *Temp = x->ObtainReadTempData("PWxsq");
		Vector *PWxsq = Temp->GetSharedElement();

		const double *xPtr = x->ObtainReadData();
		const double *tmp = nullptr;
		const double *PWxsqPtr = PWxsq->ObtainReadData();
		double *egfPtr = egf->ObtainWriteEntireData();
		integer length = PWxsq->Getlength();
		for (integer i = 0; i < egf->Getlength(); i++)
		{
			tmp = W[i]->ObtainReadData();
			egfPtr[i] = 4 * xPtr[i] * ddot_(&length, const_cast<double *> (tmp), &GLOBAL::IONE, const_cast<double *> (PWxsqPtr), &GLOBAL::IONE);
		}
	};

	void SphereConvexHull::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const SharedSpace *Temp = x->ObtainReadTempData("PWxsq");
		Vector *PWxsq = Temp->GetSharedElement();
		const double *PWxsqPtr = PWxsq->ObtainReadData();
		const double *etaxptr = etax->ObtainReadData();
		double *exixptr = exix->ObtainWriteEntireData();
		const double *tmp = nullptr;

		integer length = PWxsq->Getlength();
		for (integer i = 0; i < exix->Getlength(); i++)
		{
			tmp = W[i]->ObtainReadData();
			exixptr[i] = 4 * etaxptr[i] * ddot_(&length, const_cast<double *> (tmp), &GLOBAL::IONE, const_cast<double *> (PWxsqPtr), &GLOBAL::IONE);
		}

		Vector *xeta = x->ConstructEmpty();
		x->CopyTo(xeta);
		double *xetaPtr = xeta->ObtainWritePartialData();
		for (integer i = 0; i < xeta->Getlength(); i++)
			xetaPtr[i] *= etaxptr[i];

		Vector *Wxeta = W[0]->ConstructEmpty();
		Mani->ScaleTimesVector(x, xetaPtr[0], W[0], Wxeta);
		for (integer i = 1; i < lengthW; i++)
		{
			Mani->scalarVectorAddVector(x, xetaPtr[i], W[i], Wxeta, Wxeta);
		}
		delete xeta;
		Vector *PWxeta = Wxeta->ConstructEmpty();
		if (Hv == NULL)//--- nullptr)
			Wxeta->CopyTo(PWxeta);
		else
			(solver->*Hv)(Wxeta, PWxeta);

		const double *PWxetaPtr = PWxeta->ObtainReadData();
		const double *xPtr = x->ObtainReadData();
		for (integer i = 0; i < exix->Getlength(); i++)
		{
			tmp = W[i]->ObtainReadData();
			exixptr[i] += 8 * xPtr[i] * ddot_(&length, const_cast<double *> (tmp), &GLOBAL::IONE, const_cast<double *> (PWxetaPtr), &GLOBAL::IONE);
		}
		delete Wxeta;
		delete PWxeta;
	};

}; /*end of ROPTLIB namespace*/
