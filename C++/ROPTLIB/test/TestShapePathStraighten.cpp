//
//  TestShapePathStraighten.cpp
//  Updated
//
//  Created by Yaqing You on 2/8/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#include "test/TestShapePathStraighten.h"
using namespace ROPTLIB;

//#define MATLAB_MEX_FILE //djkfksldj


#if !defined(MATLAB_MEX_FILE) && defined(TESTSHAPEPATHSTRAIGHTEN)

int main(void)
{
    testShapePathStraighten();
#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
    return 0;
}

void testShapePathStraighten()
{
    //generate q_1 q_2
    integer numP = 101, dim = 2, numC = 11;
    //double *C1 = new double[numP*dim], *C2 = new double[numP*dim];
    double *q1 = new double[numP*dim], *q2 = new double[numP*dim];
    
    //****************************READ C1 C2*******************************
    // C1 C2 closed
    std::ifstream infile("Data_new.txt");
    if (! infile) {
        std::cout << "File did not open"<<std::endl;
    }
    else
    {
        for (integer j = 0; j < dim; j++) {
            for (integer i = 0; i < numP; i++) {
//                if (i == numP-1) {
//                    C1[i+j*numP] = C1[j*numP];
//                }
//                else
//                {
                    infile >> q1[i+j*numP];
//                }
            }
        }
        for (integer j = 0; j < dim; j++) {
            for (integer i = 0; i < numP; i++) {
//                if (i == numP-1) {
//                    C2[i+j*numP] = C2[j*numP];
//                }
//                else
//                {
                    infile >> q2[i+j*numP];
//                }
            }
        }
    }
    

// *************************DEBUG*************************
//    ForDebug::Print("C1", C1, numP*dim);
//    ForDebug::Print("q1", q1, numP*dim);
//    std::cout<<"hdfhskdjhfkjasdhfkjsdhfkahjsdfkahsdkfhasdkfhakshfaks"<<std::endl;
//    ForDebug::Print("q2", q2, numP*dim);
//    ForDebug::Print("C2", C2, numP*dim);
// *************************DEBUG*************************

    
    
    //****************************CHANGE TO q1 q2****************************
    //CurveToQ(C1, dim, numP, q1, 1);
    //CurveToQ(C2, dim, numP, q2, 1);
    
    
    integer numofmanis = 3;
    integer numofmani1 = 1;
    integer numofmani2 = 1;
    integer numofmani3 = 1;
    
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
    // ######### Read l ###########
    std::ifstream infile_2("Data_l.txt");
    if (! infile_2)
    {
        std::cout << "File did not open"<<std::endl;
    }
    else
    {
        for (integer i = 0; i < numP; i++)
        {
            infile_2 >> lOm[i];
        }
    }
    // ######### Read O ###########
    std::ifstream infile_3("Data_O.txt");
    if (! infile_3)
    {
        std::cout << "File did not open"<<std::endl;
    }
    else
    {
        for (integer i = 0; i < dim; i++)
        {
            for (integer j = 0; j < dim; j++)
            {
                infile_3 >> lOm[numP+j*dim+i];
            }
        }
    }
    // ######### Read m ###########
    std::ifstream infile_4("Data_m.txt");
    if (! infile_4)
    {
        std::cout << "File did not open"<<std::endl;
    }
    else
    {
            infile_4 >> lOm[numP+dim*dim];
    }

    
    ShapePathStraighten Prob(q1, q2, numP, dim, numC);
    Prob.SetDomain(Domain);
    
//    Prob.f(Xinitial);
    
//    RSD *RSDsolver = new RSD(&Prob, Xinitial);
//    RSDsolver->Debug = ITERRESULT;
//    RSDsolver->Max_Iteration = 50;
//    // RSDsolver->Stop_Criterion =   //user manual
//    RSDsolver->CheckParams();
//    RSDsolver->Run();
//    delete RSDsolver;
    
    LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, Xinitial);
    LRBFGSsolver->Debug = ITERRESULT;
    LRBFGSsolver->Max_Iteration = 5;
    LRBFGSsolver->CheckParams();
    LRBFGSsolver->Run();
    delete LRBFGSsolver;
    delete [] q1;
    delete [] q2;
    
    delete ll;
    delete OO;
    delete mm;
    
    delete llv;
    delete OOv;
    delete mmv;
    
    delete Domain;
    delete Xinitial;
    
}

#endif

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    printf("started===========\n");
    if (nrhs < 5) {
        mexErrMsgTxt("The number of arguments should be at least 5.\n");
    }
    
    double *C1, *C2;
    integer numP, dim, numC;
    
    C1 = mxGetPr(prhs[0]);
    C2 = mxGetPr(prhs[1]);
    numP = static_cast<integer>(mxGetScalar(prhs[2]));
    dim = static_cast<integer>(mxGetScalar(prhs[3]));
    numC = static_cast<integer>(mxGetScalar(prhs[4]));

//===================Debug========================
//    ForDebug::Print("C1", C1, numP*dim);
//    ForDebug::Print("C2", C2, numP*dim);
//    printf("\n");
//    printf("numP %d dim %d numC %d", numP, dim, numC);
//    printf("\n");
//===================Debug========================
    
    integer numofmanis = 3;
    integer numofmani1 = 1;
    integer numofmani2 = 1;
    integer numofmani3 = 1;
    
    L2Sphere *ll = new L2Sphere(numP);
    OrthGroup *OO = new OrthGroup(dim);
    Euclidean *mm = new Euclidean(1);
    
    L2SphereVariable *llv = new L2SphereVariable(numP);
    OrthGroupVariable *OOv = new OrthGroupVariable(dim);
    EucVariable *mmv = new EucVariable(1);
    
    ProductManifold *Domain = nullptr;
    Domain = new ProductManifold(numofmanis, ll, numofmani1, OO, numofmani2, mm, numofmani3);
    
    ProductElement *Xinitial = new ProductElement(numofmanis, llv, numofmani1, OOv, numofmani2, mmv, numofmani3);
    
//============Set Parameter to Get Initial Oml=============
    double w = 0.01;
    integer rotated = 1, isclosed = 1, onlyDP = 0, skipm = 4;
    char methodname[] = "LRBFGS";
    integer autoselectC = 0;
    bool swap;
    double *fopts = new double[5];
    double *comtime = new double[5];
    integer ns, lms;
    
    DriverElasticCurvesRO(C1, C2, dim, numP, w, rotated, isclosed, onlyDP, skipm, methodname, autoselectC, Xinitial, swap, fopts, comtime, ns, lms);
    
//    Xinitial->Print("Xinitial");
    printf("\n");
    printf("pointerpointer\n");
//========================================================================
    double *Xinitialptr = Xinitial->ObtainWritePartialData();
    double *gamma = new double[numP];
    dcopy_(&numP, const_cast<double *>(Xinitialptr), &GLOBAL::IONE, gamma, &GLOBAL::IONE);
//    double *gamma_coefs = new double[4*(numP-1)];
//    double *Dgamma_coefs = new double[3*(numP-1)];
    double *gamma_p = new double[numP];
//    Spline::SplineUniformPeriodic(gamma, numP, 1.0/(numP-1), gamma_coefs);
//    Spline::FirstDeri(gamma_coefs, numP, Dgamma_coefs);
//    for (integer i = 0; i < numP; i++)
//    {
//        gamma_p[i] = Spline::ValFirstDeriUniform(Dgamma_coefs, numP, 1.0/(numP-1), gamma[i]);
//    }
//    for (integer i = 0; i < numP; i++)
//    {
//        gamma_p[i] = sqrt(gamma_p[i]);
//    }
    if (isclosed)
        GradientPeriod(gamma, numP, 1.0/(numP-1), gamma_p);
    else
        Gradient(gamma, numP, 1.0/(numP-1), gamma_p);
    for (integer j = 0; j < numP; j++) {
        gamma_p[j] = sqrt(gamma_p[j]);
    }
//========================================================================
//    ForDebug::Print("Xinitial22222", Xinitialptr, numP+dim*dim+1);
//    printf("swap %d ns %d lms %d", swap, ns, lms);

//============Set Parameter to Get Initial Oml=============
    
//========consider if swap is true or false==========

    double *gamma_inverse = new double[numP];
    double *QNO = new double[dim*dim];
    double *QNm = new double[1];
    integer DIM2 = dim*dim;
    double tempdouble;
    dcopy_(&DIM2, const_cast<double *>(Xinitialptr+numP), &GLOBAL::IONE, QNO, &GLOBAL::IONE);
    dcopy_(&GLOBAL::IONE, const_cast<double *>(Xinitialptr+numP+dim*dim), &GLOBAL::IONE, QNm, &GLOBAL::IONE);
    if (swap == true){
        GammaInverse(gamma_p, numP, gamma_inverse);
        QNm[0] = 1.0 - QNm[0];
        for (integer i = 0; i < dim; i++) {
            for (integer j = i; j < dim; j++) {
                tempdouble = QNO[j*dim+i];
                QNO[j*dim+i] = QNO[i*dim+j];
                QNO[i*dim+j] = tempdouble;
            }
        }
    }
    else {
        dcopy_(&numP, const_cast<double *>(gamma_p), &GLOBAL::IONE, gamma_inverse, &GLOBAL::IONE);
    }
    dcopy_(&numP, const_cast<double *>(gamma_inverse), &GLOBAL::IONE, Xinitialptr, &GLOBAL::IONE);
    dcopy_(&DIM2, const_cast<double *>(QNO), &GLOBAL::IONE, Xinitialptr+numP, &GLOBAL::IONE);
    dcopy_(&GLOBAL::IONE, const_cast<double *>(QNm), &GLOBAL::IONE, Xinitialptr+DIM2, &GLOBAL::IONE);
    
//========consider if swap is true or false==========
    
//===========turn C1 C2 to q1 q2==============
    double *q1 = new double[numP*dim], *q2 = new double[numP*dim];
    CurveToQ(C1, dim, numP, q1, isclosed);
    CurveToQ(C2, dim, numP, q2, isclosed);
//    ForDebug::Print("q1", q1, numP*dim);
//    ForDebug::Print("q2", q2, numP*dim);
//===========turn C1 C2 to q1 q2==============

    ShapePathStraighten Prob(q1, q2, numP, dim, numC);
    Prob.SetDomain(Domain);
    
    LRBFGS *LRBFGSsolver = new LRBFGS(&Prob, Xinitial);
    LRBFGSsolver->Debug = ITERRESULT;
    LRBFGSsolver->Max_Iteration = 100;
    LRBFGSsolver->CheckParams();
    LRBFGSsolver->Run();
    
    //Flag Flag
//     double *FinalDis = (double *) mxGetData(plhs[1]);
    double FinalDis;
    FinalDis = LRBFGSsolver->Getfinalfun();
    
    
    
    
    Prob.finalPSCV->ObtainWritePartialData();
//     ForDebug::Print("log_map", Prob.log_map, numP, dim);
    mexProblem::ObtainMxArrayFromElement(plhs[0], Prob.finalPSCV);
    //plhs[0] is [2222x1 double],要不要改成矩阵形式？matlab reshape
    plhs[1] = mxCreateDoubleScalar(FinalDis);
//    plhs[1] = mxCreateDoubleScalar(static_cast<double> (FinalDis));
//    plhs[1] = FinalDis;
    
    delete LRBFGSsolver;
    delete Domain;
    delete Xinitial;
    delete [] q1;
    delete [] q2;
    delete [] fopts;
    delete [] comtime;
    delete [] gamma;
    delete [] gamma_p;
    delete [] gamma_inverse;
    delete [] QNO;
    delete [] QNm;
    /////////////////////
    
    delete ll;
    delete OO;
    delete mm;
    delete llv;
    delete OOv;
    delete mmv;
    
    return;
}

#endif
