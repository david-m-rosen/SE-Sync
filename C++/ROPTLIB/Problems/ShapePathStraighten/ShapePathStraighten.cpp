//
//  ShapePathStraighten.cpp
//  headertry
//
//  Created by Yaqing You on 10/28/15.
//  Copyright © 2015 Yaqing You. All rights reserved.
//

#include "Problems/ShapePathStraighten/ShapePathStraighten.h"

/*Define the namespace*/
namespace ROPTLIB{

	//Make sure q1 and q2 are closed
	ShapePathStraighten::ShapePathStraighten(double *inq1, double *inq2, integer innumP, integer indim, integer innumC)
	{
		q1 = inq1;
		numP = innumP;
		dim = indim;
		numC = innumC;
		q2_coefs = new double[4*dim*(numP-1) + 3*dim*(numP-1)];
		dq2_coefs = q2_coefs + 4*dim*(numP-1);
		for (integer i = 0; i < dim; i++)
		{
			Spline::SplineUniformPeriodic(inq2+i*numP, numP, 1.0/(numP-1), q2_coefs+i*4*(numP-1));
		}
		for (integer i = 0; i < dim; i++)
		{
			Spline::FirstDeri(q2_coefs+i*4*(numP-1), numP, dq2_coefs+i*3*(numP-1));
		}
        //===================== For Splines =====================
        finalPSCV = nullptr;
        log_map = new double[numP*dim];
        //===================== For Splines =====================
	};

	ShapePathStraighten::~ShapePathStraighten(void)
	{
        if(finalPSCV != nullptr)
        delete finalPSCV;
        delete [] log_map;
        delete [] q2_coefs;
	};

	double ShapePathStraiLinesearchInput(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *prob, const Solvers *solver)
	{
		return 1;
	}

	double ShapePathStraighten::f(Variable *x) const
	{
		const double *l = x->ObtainReadData();
		const double *O = l + numP;
		const double *m = O + dim * dim;
		//double *q2_new = new double[numP*dim];
        SharedSpace *Shared_q2_new = new SharedSpace(2, numP, dim); //===================
        double *q2_new = Shared_q2_new->ObtainWriteEntireData();   //==================

    
        Apply_Oml(O, m, l, numP, dim, q2_coefs, q2_new);
		//ForDebug::Print("q2q2=============", q2_new, numP, dim);
    
		//******************************compute path_x*********************************
		//initiate
		PSCVariable PSCV(numP, dim, numC);
		PSCV.Generate(q1, q2_new);
		//PSCV.Print();
		PreShapeCurves PSCurves(numP, dim, numC);
		PreShapePathStraighten PreSPSprob(numP, dim, numC);
		PreSPSprob.SetDomain(&PSCurves);
    
		//get updated preshape geodesic
		RSD *RSDsolver = new RSD(&PreSPSprob, &PSCV);
		//RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
		RSDsolver->Debug= NOOUTPUT;
		RSDsolver->LineSearch_LS = INPUTFUN;
		RSDsolver->LinesearchInput = &ShapePathStraiLinesearchInput;
		RSDsolver->Max_Iteration = 16;
        RSDsolver->IsPureLSInput = true;
		//RSDsolver->Stop_Criterion = GRAD_F;
	//    RSDsolver->CheckParams();
		RSDsolver->Run();
    
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% store the geodesic %%%%%%%%%%%%%%%%%%%%%%
        //===================== For Splines =====================
        if(finalPSCV == nullptr)
            finalPSCV = PSCV.ConstructEmpty();
        RSDsolver->GetXopt()->CopyTo(finalPSCV);
        //===================== For Splines =====================
        
		//store preshape geodesic to compute distance and eta
		//current preshape geodesic between q1 and q2_new stored in PreGeodesic
		PSCVariable *PreGeodesic = PSCV.ConstructEmpty();
		RSDsolver->GetXopt()->CopyTo(PreGeodesic);
    
		//PreGeodesic->Print();
    
    
		//#############################Compute Distance##############################
		  //######Compute Dalpha######
		const double *GeoPath = PreGeodesic->ObtainReadData();
		//PreGeodesic->Print("kdsjfl;akjsdlfjals;dfjal;dsjfal;dsjfl;asdjfl;asdjfl;asjd");
		double *Dalpha = new double[numP*dim*numC];
		integer stp = numC - 1;
		for (integer t = 0; t < numC; t++)
		{
			if (t != 0) {
				for (integer j = 0; j < dim; j++)
				{
					for (integer i = 0; i < numP; i++)
					{
						Dalpha[t*numP*dim+j*numP+i] = stp*(GeoPath[t*numP*dim+j*numP+i] - GeoPath[(t-1)*numP*dim+j*numP+i]);
					}
				}
			}
			//Project c(tau/numC) into T_alpha(M)
		   if (t == 0)
		   {
			   for (integer j = 0; j < dim; j++)
			   {
					for (integer i = 0; i < numP; i++)
					{
						Dalpha[j*numP+i] = 0.0;
					}
				}
			}
			else
			{
				PreShapePathStraighten::Item_2(GeoPath + t*numP*dim, numP, dim, Dalpha + t*numP*dim);
			}
		}
		//ForDebug::Print("alphaalphalphalphalphalph", Dalpha, numP, dim, numC);
    
		//#######Compute Distance#########
		double intv;
		double *temp = new double[numC];
		double distance;
    
		for (integer i = 0; i < numC; i++)
		{
			temp[i] = PreShapePathStraighten::InnerProd_Q(Dalpha+i*numP*dim, Dalpha+i*numP*dim, numP, dim);
			temp[i] = sqrt(temp[i]);
		}
		intv = 1.0/(numC-1);
		distance = ElasticCurvesRO::Trapz(temp, numC, intv);
		//x->Print("================================");
    
		//#############################Compute and store eta##############################
		integer NXD = numP*dim;
		double coeff;
		SharedSpace *SharedEta = new SharedSpace(2, numP, dim);
		double *eta = SharedEta->ObtainWriteEntireData();
		dcopy_(&NXD, const_cast<double *>(Dalpha+(numC-1)*numP*dim), &GLOBAL::IONE, eta, &GLOBAL::IONE);
		coeff = 1.0/std::sqrt(PreShapePathStraighten::InnerProd_Q(Dalpha+(numC-1)*numP*dim, Dalpha+(numC-1)*numP*dim, numP, dim));
		dscal_(&NXD, &coeff, eta, &GLOBAL::IONE);
    
        //===================== For Splines =====================
        //===================== For Splines =====================
//        if(log_map == nullptr)
//            log_map = new double[numP*dim];
        dcopy_(&NXD, const_cast<double *>(Dalpha+numP*dim), &GLOBAL::IONE, log_map, &GLOBAL::IONE);
        
        
        

		x->AddToTempData("eta", SharedEta);
        x->AddToTempData("q2_new", Shared_q2_new);
    
		//############################# Return ##############################
        

		delete RSDsolver;
		//delete [] q2_new;
		delete [] Dalpha;
		delete [] temp;
        delete PreGeodesic;
		return distance;
	}



	void ShapePathStraighten::EucGrad(Variable *x, Vector *egf) const
	{
		const double *eta = x->ObtainReadTempData("eta")->ObtainReadData();
		double *egfGrad = egf->ObtainWriteEntireData();
    
    
		//const SharedSpace *Temp = x->ObtainReadTempData("Dalpha");
		//Temp->GetSharedElement()->CopyTo(egf);
		const double *l = x->ObtainReadData();
		const double *O = l + numP;
		const double *m = O + dim*dim;
		double *l_temp = new double[numP];
		double *gamma = new double[numP];
		double intv = 1.0/(numP-1);
		integer NXD = numP*dim;
		double *q2_new = new double[NXD];
		double *q2_temp = new double[NXD];
		integer innumP = numP;
		integer indim = dim;
		// integer innumC = numC;
    
		//Compute int_0^1 l^2 + m mod 1
		for (integer i = 0; i < numP; i++)
		{
			l_temp[i] = l[i]*l[i];
		}
		ElasticCurvesRO::CumTrapz(l_temp, numP, intv, gamma);
		//Get m0
		double m0 = std::fmod(m[0], 1);
		while (m0 > 1)
		{
			m0 -= 1.0;
		}
		while (m0 < 0)
		{
			m0 += 1.0;
		}
		//#### get updated \gamma ####
		for (integer i = 0; i < numP; i++)
		{
			gamma[i] += m0;
			gamma[i] = (gamma[i] > 1) ? gamma[i] - 1 : gamma[i];
		}
		//#### get q2 \circ \gamma ####
		for (integer i = 0; i < numP; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				q2_new[i+j*numP] = Spline::ValSplineUniform(q2_coefs+j*4*(numP-1), numP, 1.0/(numP-1), gamma[i]);
			}
		}
    
		//#### l * q2_new ####
		double coeff_temp;
		for (integer i = 0; i < numP; i++)
		{
			coeff_temp = fabs(l[i]);
			for (integer j = 0; j < dim; j++)
			{
				q2_new[j*numP+i] = q2_new[j*numP+i] * coeff_temp;
			}
		}

		//################## A = int(eta * q2_new) ##################
		double *A = new double[dim*dim];
		double *product = new double[numP];
		for (integer i = 0; i < dim; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				//### eta_i * q_j ###
				for (integer k = 0; k < numP; k++)
				{
					product[k] = eta[k+i*numP] * q2_new[k+j*numP];
				}
				//ForDebug::Print("product", product, numP);
				A[j*dim+i] = ElasticCurvesRO::Trapz(product, numP, intv);
			}
		}
    
		//ForDebug::Print("A", A, dim, dim);
    
    
		//################## calculate x(as xx, to avoid duplicates) ##################
		for (integer i = 0; i < numP; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				q2_new[i+j*numP] = Spline::ValSplineUniform(q2_coefs+j*4*(numP-1), numP, 1.0/(numP-1), gamma[i]);
			}
		}
    
		// result <- q2_new * O^T
		//details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		char *transn = const_cast<char *> ("n");
		char *transt = const_cast<char *> ("t");
		for (integer i = 0; i < NXD; i++)
		{
			q2_temp[i] = 0.0;
		}
		dgemm_(transn, transt, &innumP, &indim, &indim, &GLOBAL::DONE, const_cast<double *>(q2_new), &innumP, const_cast<double *>(O), &indim, &GLOBAL::DZERO, q2_temp, &innumP);
    
		//############### calculate <eta , q2_temp>
		double *xx = new double[numP];
		//initial xx
		for (integer i = 0; i < numP; i++)
		{
			xx[i] = 0.0;
		}
		for (integer i = 0; i < numP; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				xx[i] += eta[j*numP+i] * q2_temp[j*numP+i];
			}
		}
    
    
		//############################# calculate y' ##############################
		double *q2p = new double[NXD];
		for (integer i = 0; i < numP; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				q2p[i+j*numP] = Spline::ValFirstDeriUniform(dq2_coefs+j*3*(numP-1), numP, 1.0/(numP-1), gamma[i]);
			}
		}

		// result <- q2_new * O^T
		//details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		for (integer i = 0; i < NXD; i++)
		{
			q2_temp[i] = 0.0;
		}
		dgemm_(transn, transt, &innumP, &indim, &indim, &GLOBAL::DONE, const_cast<double *>(q2p), &innumP, const_cast<double *>(O), &indim, &GLOBAL::DONE, q2_temp, &innumP);
    

    
		for (integer i = 0; i < numP; i++)
		{
			coeff_temp = fabs(l[i]);
			for (integer j = 0; j < dim; j++)
			{
				q2p[j*numP+i] = q2_temp[j*numP+i] * coeff_temp;
			}
		}
    
		//############### calculate <eta , q2p>
		double *yp = new double[numP];
		for (integer i = 0; i < numP; i++)
		{
			yp[i] = 0.0;
		}
		for (integer i = 0; i < numP; i++)
		{
			for (integer j = 0; j < dim; j++)
			{
				yp[i] += eta[j*numP+i] * q2p[j*numP+i];
			}
		}
    
		//######################### GRAD_m ##########################
		double Grad_m;
		Grad_m = ElasticCurvesRO::Trapz(yp, numP, intv);
    
		//######################### Grad_l = xx-2yl ##########################
		double *y = new double[numP];
		double *Grad_l = new double[numP];
		ElasticCurvesRO::CumTrapz(yp, numP, intv, y);
		for (integer i = 0; i < numP; i++)
		{
			Grad_l[i] = xx[i] - 2*y[i]*l[i];
		}
    
		//#################### Grad_O = A #####################
    
		//#################### Copy Grads to egfGrad #####################
		for (integer i = 0; i < numP; i++)
		{
			egfGrad[i] = Grad_l[i];
		}
		for (integer i = 0; i < dim*dim; i++)
		{
			egfGrad[numP+i] = A[i];
		}
		egfGrad[numP+dim*dim] = Grad_m;
    
		//######### Delete #########
		delete [] l_temp;
		delete [] gamma;
		delete [] q2_new;
		delete [] q2_temp;
		delete [] A;
		delete [] product;
		delete [] xx;
		delete [] q2p;
		delete [] yp;
		delete [] y;
		delete [] Grad_l;
    
    
	};



	void ShapePathStraighten::Apply_Oml(const double *O, const double *m, const double *l, integer innumP, integer indim, const double *q2_coeff, double *q2_new)
	{
		char *transn = const_cast<char *> ("n");
		char *transt = const_cast<char *> ("t");
		double *gamma_temp = new double[innumP];
		int NXD = innumP*indim;
		double intv = 1.0/(innumP-1);
		double *l_temp = new double[innumP];
		double *q2_temp = new double[NXD];
		double coeff_temp;
    
		//Compute l.^2
		for (integer i = 0; i < innumP; i++)
		{
			l_temp[i] = l[i]*l[i];
		}
    
		ElasticCurvesRO::CumTrapz(l_temp, innumP, intv, gamma_temp);
		//get m0
		double m0 = std::fmod(m[0], 1);
		while (m0 > 1)
		{
			m0 -= 1.0;
		}
		while (m0 < 0)
		{
			m0 += 1.0;
		}
		//int_(l^2) mod 1
		for (integer i = 0; i < innumP; i++)  //mod 1
		{
			gamma_temp[i] += m0;
			gamma_temp[i] = (gamma_temp[i] > 1) ? gamma_temp[i] - 1 : gamma_temp[i];
		}
    
		//ppval
		for (integer i = 0; i < innumP; i++)
		{
			for (integer j = 0; j < indim; j++)
			{
				q2_new[i+j*innumP] = Spline::ValSplineUniform(q2_coeff+j*4*(innumP-1), innumP, 1.0/(innumP-1), gamma_temp[i]);
			}
		}
		// result <- q2_new * O^T
		//details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
    
		for (integer i = 0; i < NXD; i++)
		{
			q2_temp[i] = 0.0;
		}
		dgemm_(transn, transt, &innumP, &indim, &indim, &GLOBAL::DONE, const_cast<double *>(q2_new), &innumP, const_cast<double *>(O), &indim, &GLOBAL::DZERO, q2_temp, &innumP);
    
		for (integer i = 0; i < innumP; i++)
		{
			coeff_temp = fabs(l[i]);
			for (integer j = 0; j < indim; j++)
			{
				q2_new[j*innumP+i] = q2_temp[j*innumP+i] * coeff_temp;
			}
		}
    
		delete [] gamma_temp;
		delete [] l_temp;
		delete [] q2_temp;
	}

}; /*end of ROPTLIB namespace*/


