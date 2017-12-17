
#include "Problems/Problem.h"

/*Define the namespace*/
namespace ROPTLIB{

	void Problem::CheckGradHessian(const Variable *xin) const
	{
		UseGrad = true;
		UseHess = true;
		integer length;
		double normxi;
		double t, fx, fy;
		double *X, *Y;
		Vector *etax;
		Variable *x = xin->ConstructEmpty();
		xin->CopyTo(x);
		if (Domain->GetIsIntrinsic())
			etax = Domain->GetEMPTYINTR()->ConstructEmpty();
		else
			etax = Domain->GetEMPTYEXTR()->ConstructEmpty();
		etax->RandUnform();
		Vector *xi = etax->ConstructEmpty();
		Vector *gfx = etax->ConstructEmpty();
		Vector *Hv = etax->ConstructEmpty();
        Variable *y = x->ConstructEmpty();
        fx = f(x);
		Grad(x, gfx);
        //gfx->Print("gfx:");//---
		gfx->CopyTo(etax);//--

		//double *etaxTV = etax->ObtainWriteEntireData();///---
		//integer nnn = etax->Getlength();
		//for (integer i = 0; i < nnn; i++)//--
		//{
		//	etaxTV[i] = sin(static_cast<double> (i) / (etax->Getlength() - 1) / 2);
		//}
		//for (integer i = 0; i < 5; i++)//---
		//	etaxTV[nnn - 1 - i] = 0;//--

		//etax->Print("etax:");//--
        Domain->Projection(x, etax, xi);
        normxi = sqrt(Domain->Metric(x, xi, xi));
		Domain->ScaleTimesVector(x, 100 / normxi, xi, xi); // initial length of xi is 100
		//xi->Print("xi:");//---
		// the length of xi variances from 100 to 100*2^(-35) approx 6e-9
		t = 1;
		length = 35;
		X = new double[length * 2]; 
		Y = X + length;
		for (integer i = 0; i < length; i++)
		{
			Domain->Retraction(x, xi, y, 0);
			fy = f(y);
			//y->Print("y:");//----
			HessianEta(x, xi, Hv);
			Y[i] = log(fabs(fy - fx - Domain->Metric(x, gfx, xi) - 0.5 * Domain->Metric(x, xi, Hv)));
			X[i] = 0.5 * log(Domain->Metric(x, xi, xi));
			printf("i:%d,|eta|:%.3e,(fy-fx)/<gfx,eta>:%.3e,(fy-fx-<gfx,eta>)/<0.5 eta, Hessian eta>:%.3e\n", i,
				sqrt(Domain->Metric(x, xi, xi)), (fy - fx) / Domain->Metric(x, gfx, xi),
				(fy - fx - Domain->Metric(x, gfx, xi)) / (0.5 * Domain->Metric(x, xi, Hv)));
			Domain->ScaleTimesVector(x, 0.5, xi, xi);
		}

		printf("CHECK GRADIENT:\n");
		printf("\tSuppose the point is not a critical point.\n");
		printf("\tIf there exists an interval of |eta| such that (fy - fx) / <gfx, eta>\n");
		printf("\tapproximates ONE, then the gradient is probably correct!\n");

		printf("CHECK THE ACTION OF THE HESSIAN (PRESUME GRADIENT IS CORRECT):\n");
		printf("\tSuppose the retraction is second order or the point is a critical point.\n");
		printf("\tIf there exists an interval of |eta| such that (fy-fx-<gfx,eta>)/<0.5 eta, Hessian eta>\n");
		printf("\tapproximates ONE, then the action of Hessian is probably correct.\n");

		////TEST IDEA2: 
		//for (integer i = 1; i < length - 1; i++)
		//	printf("log(|eta|):%.3e, slope:%.3e\n", X[i], (Y[i + 1] - Y[i - 1]) / (X[i + 1] - X[i - 1]));
		//printf("CHECK GRADIENT:\n");
		//printf("\tIf there exists an interval of |eta| such that the slopes \n");
		//printf("\tapproximate TWO, then the gradient is probably correct!\n");

		//printf("CHECK THE ACTION OF THE HESSIAN (PRESUME GRADIENT IS CORRECT AND\n");
		//printf("THE COST FUNCTION IS NOT ONLY QUADRATIC):\n");
		//printf("\tIf there exists an interval of |eta| such that the slopes\n");
		//printf("\tapproximate THREE, then the action of Hessian is probably correct.\n");

		//x->Print("1, x:", false);//---
		delete xi;
		//x->Print("2, x:", false);//---
		//gfx->Print("2, gfx:", false);//---
		delete gfx;
		//x->Print("3, x:", false);//---
		delete y;
		//x->Print("4, x:", false);//---
		delete Hv;
		//x->Print("5, x:", false);//---
		delete[] X;
		delete etax;
		//x->Print("x:", false);//---
		delete x;
	};

	void Problem::Grad(Variable *x, Vector *gf) const
	{
		if (!Domain->GetIsIntrinsic())
		{
			RieGrad(x, gf);
			return;
		}
		Vector *exgf = Domain->GetEMPTYEXTR()->ConstructEmpty();
		RieGrad(x, exgf);
		//exgf->Print("exgf:");//---
		Domain->ObtainIntr(x, exgf, gf);
		delete exgf;
	};

	void Problem::HessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		if (!Domain->GetIsIntrinsic())
		{
			RieHessianEta(x, etax, xix);
			return;
		}

		Vector *exxix = Domain->GetEMPTYEXTR()->ConstructEmpty();
		Vector *exetax = Domain->GetEMPTYEXTR()->ConstructEmpty();
		Domain->ObtainExtr(x, etax, exetax);
		RieHessianEta(x, exetax, exxix);
		Domain->ObtainIntr(x, exxix, xix);
		delete exxix;
		delete exetax;
	}

	void Problem::PreConditioner(Element *x, Element *inVec, Element *outVec) const
	{
		// default one means no preconditioner.
		inVec->CopyTo(outVec);
	};

	void Problem::RieGrad(Variable *x, Vector *gf) const
	{
		EucGrad(x, gf);		

		///*For some of the manifolds, converting the Euclidean gradient to intrinsic representation
		//is equivalent to converting the Euclidean gradeitn to Riemannian gradient and then converting
		//the Riemannian gradient to the intrinsic representation. Therefore, one can avoid computing
		//the Riemannian gradient to save computations.*/
		//if (!Domain->GetIsIntrinsic() || UseHess || 
		//	Domain->GetName() == "SPDTensor" || Domain->GetName() == "SPDManifold" ||
		//	Domain->GetName() == "PreShapeCurves" || Domain->GetName() == "LowRank" ||
		//	Domain->GetName() == "EucPositive" || Domain->GetName() == "ElasticShape" ||
		//	Domain->GetName() == "CpxNStQOrth")
		//{
		//	Domain->EucGradToGrad(x, gf, gf, this);
		//}

		/*using this function for simplicity. Above comparison statement is too long and is inefficient for
		small size problems.*/
		Domain->EucGradToGrad(x, gf, gf, this);
	};

	void Problem::RieHessianEta(Variable *x, Vector *etax, Vector *xix) const
	{
		EucHessianEta(x, etax, xix);
		Domain->EucHvToHv(x, etax, xix, xix, this);
	};

	void Problem::EucGrad(Variable *x, Vector *egf) const
	{
		printf("Euclidean Gradient has not been done!\n");
	};

	void Problem::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		/*
		finite difference to approximate the action of the Hessian.
		Since everything is done in Euclidean space, no retraction and vector transport are necessary.
		*/
		Variable *y = x->ConstructEmpty();
		Vector *gfy = etax->ConstructEmpty();
		double normetax = sqrt(Domain->Metric(x, etax, etax));
		double factor = 1e-5 / normetax;
		Domain->ScaleTimesVector(x, factor, etax, exix);
		Domain->VectorAddVector(x, x, exix, y);

		/*
		f(y) is evaluated before EucGrad since some computations, which are needed in EucGrad, are done in f(y).
		The Euclidean gradient uses extrinsic approach.
		*/
		f(y); EucGrad(y, gfy);
		const SharedSpace *Sharedegf = x->ObtainReadTempData("EGrad");
		Vector *gfx = Sharedegf->GetSharedElement();
		Domain->VectorLinearCombination(x, 1.0 / factor, gfy, -1.0 / factor, gfx, exix); //

		delete y;
		delete gfy;
	};

	void Problem::SetDomain(Manifold *inDomain)
	{
		Domain = inDomain;
	};

	Problem::~Problem(void)
	{
	};

	void Problem::SetUseGrad(bool usegrad) const
	{
		UseGrad = usegrad;
	};

	void Problem::SetUseHess(bool usehess) const
	{
		UseHess = usehess;
	};
}; /*end of ROPTLIB namespace*/
