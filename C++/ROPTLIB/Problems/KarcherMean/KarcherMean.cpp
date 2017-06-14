//
//  KarcherMean.cpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#include "Problems/KarcherMean/KarcherMean.h"

/*Define the namespace*/
namespace ROPTLIB{
    KarcherMean::KarcherMean(double *Curves, integer innumP, integer indim, integer innumC, integer innumS)
    {
        Qs = Curves;
        numP = innumP;
        dim = indim;
        numC = innumC;
        numS = innumS;
    }

    KarcherMean::~KarcherMean(void)
    {
    };
    
    double LinesearchInput(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *Prob, const Solvers *solver)
    {
        return 1;
    }

    void KarcherMean::ComputeGeodesic(const double *q1, const double *q2, double &distance, double *w) const
    {
        integer numofmanis = 3;    //what is this ?
        integer numofmani1 = 1;    //what is this ?
        integer numofmani2 = 1;    //what is this ?
        integer numofmani3 = 1;    //what is this ?
        
        L2Sphere *ll = new L2Sphere(numP);
        OrthGroup *OO = new OrthGroup(dim);
        Euclidean *mm = new Euclidean(1);
        
        L2SphereVariable *llv = new L2SphereVariable(numP);
        OrthGroupVariable *OOv = new OrthGroupVariable(dim);
        EucVariable *mmv = new EucVariable(1);
        
        ProductManifold *Domain = nullptr;
        Domain = new ProductManifold(numofmanis, ll, numofmani1, OO, numofmani2, mm, numofmani3);
        
        ProductElement *Xinitial = new ProductElement(numofmanis, llv, numofmani1, OOv, numofmani2, mmv, numofmani3);
        
        //Xinitial->RandInManifold();
        double *lOm = Xinitial->ObtainWriteEntireData();
        for (integer i = 0; i < numP; i++) {
            lOm[i] = 1.0;
        }
        for (integer i = 0; i < dim; i++) {
            for (integer j = 0; j < dim; j++) {
                if(i == j) lOm[numP+j+i*dim] = 1.0;
                else lOm[numP+j+i*dim] = 0.0;
            }
        }
        lOm[numP+dim*dim] = 0.0;
        
        
        //double *lOm = Xinitial->ObtainWriteEntireData();
        // ######### Read l ###########
        // ######### Read O ###########
        // ######### Read m ###########
        
        ShapePathStraighten Prob(const_cast<double *>(q1), const_cast<double *>(q2), numP, dim, numC);
        Prob.SetDomain(Domain);
        LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, Xinitial);
        LRBFGSsolver->Debug = NOOUTPUT;
        LRBFGSsolver->Max_Iteration = 20;
        LRBFGSsolver->Stop_Criterion = FUN_REL;
        LRBFGSsolver->Tolerance = 1e-3;
        //LRBFGSsolver->CheckParams();
        LRBFGSsolver->Run();
        
        const double *q2_new = LRBFGSsolver->GetXopt()->ObtainReadTempData("q2_new")->ObtainReadData();
        
        PSCVariable PSCV(numP, dim, numC);
        PSCV.Generate(const_cast<double *>(q1), const_cast<double *>(q2_new));
        //PSCV.Print();
        //ForDebug::Print("3rd", PSCV.ObtainReadData() + 3 * numP * dim, numP, dim); //----D
        PreShapeCurves PSCurves(numP, dim, numC);
        PreShapePathStraighten PreSPSprob(numP, dim, numC);
        PreSPSprob.SetDomain(&PSCurves);
        
         RSD *RSDsolver = new RSD(&PreSPSprob, &PSCV);
         //RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
         RSDsolver->Debug= NOOUTPUT;
         RSDsolver->LineSearch_LS = INPUTFUN;
         RSDsolver->LinesearchInput = &LinesearchInput;
         RSDsolver->Max_Iteration = 10;
         RSDsolver->Stop_Criterion = GRAD_F;
         RSDsolver->Tolerance = 1e-10;
         //RSDsolver->CheckParams();
         RSDsolver->Run();
         const double *Dalpha = RSDsolver->GetXopt()->ObtainReadTempData("Dalpha")->GetSharedElement()->ObtainReadData();
         //RSDsolver->GetXopt()->Print("OPT:", false);//---
         integer NXD = numP*dim;
         dcopy_(&NXD, const_cast<double *> (Dalpha + NXD*(numC-1)), &GLOBAL::IONE, w, &GLOBAL::IONE); //======
         double coeff = -1.0;//--
         //double coeff = -1.0/std::sqrt(PreShapePathStraighten::InnerProd_Q(w, w, numP, dim));
         dscal_(&NXD, &coeff, w, &GLOBAL::IONE);
        
        
       //const double *eta = LRBFGSsolver->GetXopt()->ObtainReadTempData("eta")->ObtainReadData();
        distance = LRBFGSsolver->Getfinalfun();
        
        
        //dcopy_(&NXD, const_cast<double *> (eta), &GLOBAL::IONE, w, &GLOBAL::IONE);
        delete RSDsolver;
        delete LRBFGSsolver;

    }

    double KarcherMean::f(Variable *x) const
    {
        const double *Curve_x = x->ObtainReadData();
        //x->Print("fhjskdhashfkjasdhfklajhsdflkahskl");

        double distance = 0.0, dis_temp;
        integer NXD = numP*dim;
        double *w = new double[NXD];
        SharedSpace *SharedGrad = new SharedSpace(2, numP, dim);
        double *grad = SharedGrad->ObtainWriteEntireData();
        for (integer i = 0; i < NXD; i++) grad[i] = 0.0;  //initialize
        for (integer i = 0; i < numS; i++)
        {
            //ComputeGeodesic(Curve_x, Qs+i*numP*dim, dis_temp, w);
            ComputeGeodesic(Qs+i*numP*dim, Curve_x, dis_temp, w);    ////==========================gaigaigai=======
            distance = distance + dis_temp*dis_temp;
            daxpy_(&NXD, &GLOBAL::DNONE, w, &GLOBAL::IONE, grad, &GLOBAL::IONE);
        }
        distance = distance/(2*numS);
        double coeff = 1.0/double(numS);
        dscal_(&NXD, &coeff, grad, &GLOBAL::IONE);
        
        x->AddToTempData("grad", SharedGrad);
        
        delete[] w;
        
        return distance;
    }

    void KarcherMean::EucGrad(Variable *x, Vector *egf) const
    {
        const double *grad = x->ObtainReadTempData("grad")->ObtainReadData();
        double *egfGrad = egf->ObtainWriteEntireData();
        integer NXD = numP*dim;
        dcopy_(&NXD, const_cast<double *> (grad), &GLOBAL::IONE, egfGrad, &GLOBAL::IONE);
    }

}; /*end of ROPTLIB namespace*/
