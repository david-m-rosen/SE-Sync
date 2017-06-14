
#include "Manifolds/PreShapeCurves/PreShapeCurves.h"

/*Define the namespace*/
namespace ROPTLIB{
	PreShapeCurves::PreShapeCurves(integer r, integer c, integer n)
	{
		numP = r;
		dim = c;
		numC = n;
		IsIntrApproach = false;//tangent vector intrinsic
		HasHHR = false;
		UpdBetaAlone = false;
		name.assign("PreShapeCurves");
		IntrinsicDim = n * r * c;
		ExtrinsicDim = n * r * c;
		EMPTYEXTR = new PSCVector(r, c, n);
		EMPTYINTR = new PSCVector(r, c, n);
	};

	PreShapeCurves::~PreShapeCurves(void)
	{
		delete EMPTYEXTR;
		delete EMPTYINTR;
	};

	void PreShapeCurves::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* result, const Problem *prob) const
	{
		exix->CopyTo(result);
	};

	void PreShapeCurves::CheckParams(void) const
	{
		Manifold::CheckParams();
		printf("%s PARAMETERS:\n", name.c_str());
		if (dim == 1 && numC == 1)
			printf("numP          :%15d\n", numP);
		else
		if (numC == 1)
		{
			printf("numP          :%15d,\t", numP);
			printf("dim           :%15d\n", dim);
		}
		else
		{
			printf("numP          :%15d,\t", numP);
			printf("dim           :%15d\n", dim);
			printf("numC          :%15d\n", numC);
		}
	};

	void PreShapeCurves::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
	{
		const double *Path_x = x->ObtainReadData();
		const double *Egrad = egf->ObtainReadData();   //**************Check!!***************
		double *Rgrad = gf->ObtainWriteEntireData();
		// get numP=innumP, dim=indim, numC=innumC
		double *u = new double[numC*numP*dim];
		double *utilde = new double[numC*numP*dim];
    
		//****************************************************************
		// [u] = CovIntegral(Dalpha, alpha)
		CovIntegral(Egrad, Path_x, numC, numP, dim, u);
    
		//****************************************************************
		// [utilde] = BackTrans(u, alpha)
		BackTrans(u, Path_x, numC, numP, dim, utilde);
    
		//****************************************************************
		// [w] = GradVec(utilde, u)
		GradVec(utilde, u, numC, numP, dim, Rgrad);    //Check!!*******************************
		//ForDebug::Print("Rgrad", Rgrad, numP,dim,numC);
    
		delete[] u;
		delete[] utilde;
	};

	double PreShapeCurves::Metric(Variable *x, Vector *etax, Vector *xix) const
	{
		//const double *Path_x = x->ObtainReadData();
		const double *vec_1 = etax->ObtainReadData();  //RemannianGrad
		const double *vec_2 = xix->ObtainReadData();
		double intv;
		double *temp = new double[numC];
		double result;

		for (integer i = 0; i < numC; i++)
		{
			temp[i] = InnerProd_Q(vec_1+i*numP*dim, vec_2+i*numP*dim, numP, dim);
		}
    
		intv = 1.0/(numC-1);
		result = 0.5*ElasticCurvesRO::Trapz(temp, numC, intv);
    
		delete[] temp;
		return result;
	};

	void PreShapeCurves::Retraction(Variable *x, Vector *etax, Variable *result) const
	{
		// x is the Path, etax is update direction, result stores the updated Path.
		const double *Path_x = x->ObtainReadData();
		const double *direc = etax->ObtainReadData();
		double *Path_new = result->ObtainWriteEntireData();
		integer NXD = numP*dim;
		double *Path_temp = new double[NXD*numC];
		for (integer t = 0; t < numC; t++) {
			for (integer i = 0; i < NXD; i++) {
				Path_temp[t*NXD+i] = Path_x[t*NXD+i] + direc[t*NXD+i];
			}
			PreShapePathStraighten::Item_1(Path_temp+t*NXD, numP, dim, Path_new+t*NXD);
		}
    
		delete[] Path_temp;
	};

	//Calculating inner prod of two points q1, q2 in shape space.
	double PreShapeCurves::InnerProd_Q(const double *q1, const double *q2, integer innumP, integer indim)
	{
		double intv, result;
		double *PInnerProd = new double[innumP];
		ElasticCurvesRO::PointwiseInnerProd(q1, q2, indim, innumP, PInnerProd);
		intv = 1.0/(innumP-1);
		result = ElasticCurvesRO::Trapz(PInnerProd, innumP, intv);
		delete [] PInnerProd;
		return result;
	}



	//****************************************************************
	// [u] = CovIntegral(Dalpha, alpha)
	void PreShapeCurves::CovIntegral(const double *Dalpha, const double *alpha, integer innumC, integer innumP, integer indim, double *u)
	{
		double *ubar = new double[innumC*innumP*indim];
		double coeff;
		integer NXD = innumP*indim;
		for (integer i = 0; i < innumP*indim; i++) {
			u[i] = 0.0;
		}
		coeff = 1.0/(innumC-1);
		for (integer t = 0; t < innumC-1; t++) {
			PreShapePathStraighten::Item_3(u+t*innumP*indim, alpha+t*innumP*indim, alpha+(t+1)*innumP*indim, innumP, indim, ubar+t*innumP*indim);
			daxpy_(&NXD, &coeff, const_cast<double *> (Dalpha)+(t+1)*innumP*indim, &GLOBAL::IONE, ubar+t*innumP*indim, &GLOBAL::IONE);
			dcopy_(&NXD, ubar+t*innumP*indim, &GLOBAL::IONE, u+(t+1)*innumP*indim, &GLOBAL::IONE);
		}
		delete [] ubar;
	}


	//****************************************************************
	// [utilde] = BackTrans(u, alpha)
	void PreShapeCurves::BackTrans(const double *u, const double *alpha, integer innumC, integer innumP, integer indim, double *utilde)
	{
		double *temp = new double[innumP*indim];
		double l, coeff, coeff_1;
		integer NXD = innumP*indim;
    
		dcopy_(&NXD, const_cast<double *>(u)+(innumC-1)*NXD, &GLOBAL::IONE, utilde+(innumC-1)*NXD, &GLOBAL::IONE);
		l = std::sqrt(InnerProd_Q(u+(innumC-1)*NXD, u+(innumC-1)*NXD, innumP, indim));
		for (integer t = innumC-2; t > -1; t--) {
			dcopy_(&NXD, utilde+(t+1)*NXD, &GLOBAL::IONE, temp, &GLOBAL::IONE);
			PreShapePathStraighten::Item_2(alpha+t*NXD, innumP, indim, temp);
			coeff = std::sqrt(InnerProd_Q(temp, temp, innumP, indim));
			if (coeff < 1e-8) {
				for (integer i = 0; i < NXD; i++) {
					utilde[t*NXD + i] = 0.0;
				}
			}
			else {
				coeff_1 = l/std::sqrt(InnerProd_Q(temp, temp, innumP, indim));
				dscal_(&NXD, &coeff_1, temp, &GLOBAL::IONE);
				dcopy_(&NXD, temp, &GLOBAL::IONE, utilde+t*NXD, &GLOBAL::IONE);
			}
		}
		delete [] temp;
	}


	//****************************************************************
	// [w] = GradVec(utilde, u)
	void PreShapeCurves::GradVec(const double *utilde, const double *u, integer innumC, integer innumP, integer indim, double *w)
	{
		//ForDebug::Print("u", u, innumP,indim,innumC);
    
		double coeff;
		integer NXD = innumP*indim;
		double *temp = new double[NXD];
    
		for (integer i = 0; i < NXD; i++) {
			w[i] = 0.0;
		}
		for (integer t = 1; t < innumC; t++) {
			coeff = static_cast<double>(-t)/static_cast<double>(innumC-1);
			dcopy_(&NXD, const_cast<double *>(u)+t*NXD, &GLOBAL::IONE, temp, &GLOBAL::IONE);
			daxpy_(&NXD, &coeff, const_cast<double *>(utilde)+t*NXD, &GLOBAL::IONE, temp, &GLOBAL::IONE);
			dcopy_(&NXD, temp, &GLOBAL::IONE, w+t*NXD, &GLOBAL::IONE);
		}
		delete [] temp;
	}

}; /*end of ROPTLIB namespace*/
