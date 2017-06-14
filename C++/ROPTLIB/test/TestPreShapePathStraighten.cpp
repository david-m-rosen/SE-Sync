//
//  TestPreShapePathStraighten.cpp
//  headertry
//
//  Created by Yaqing You on 11/3/15.
//  Copyright © 2015 Yaqing You. All rights reserved.
//

#include "test/TestPreShapePathStraighten.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTPRESHAPEPATHSTRAIGHTEN)


int main(void)
{
    testPreShapePathStraighten();
    return 0;
}

double LinesearchInput(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *prob, const Solvers *solver)
{
    return 1;
}

void testPreShapePathStraighten()
{
    //generate q_1 q_2
    integer numP = 101, dim = 2, numC = 11;
    //double *C1 = new double[numP*dim], *C2 = new double[numP*dim];
    double *q1 = new double[numP*dim], *q2 = new double[numP*dim];
    

//****************************READ C1 C2*******************************
    std::ifstream infile("D:/Data.txt");
    if (! infile) {
        printf("File did not open\n");
		return;

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
    
    //ForDebug::Print("q1", q1, numP, dim);   //----D
    //ForDebug::Print("q2", q2, numP, dim);   //----D

//****************************CHANGE TO q1 q2****************************
    //CurveToQ(C1, dim, numP, q1, 1);
    //CurveToQ(C2, dim, numP, q2, 1);
    
    //ForDebug::Print("q1", q1, numP, dim); //----D
    
//*****************************COMPUTE PATH******************************
    PSCVariable PSCV(numP, dim, numC);
    PSCV.Generate(q1, q2);
    //PSCV.Print();
    //ForDebug::Print("3rd", PSCV.ObtainReadData() + 3 * numP * dim, numP, dim); //----D
    PreShapeCurves PSCurves(numP, dim, numC);
    PreShapePathStraighten PreSPSprob(numP, dim, numC);
    PreSPSprob.SetDomain(&PSCurves);
    

     RSD *RSDsolver = new RSD(&PreSPSprob, &PSCV);
     //RSDsolver->LineSearch_LS = static_cast<LSAlgo> (i);
     RSDsolver->Debug= ITERRESULT;
     RSDsolver->LineSearch_LS = INPUTFUN;
     RSDsolver->LinesearchInput = &LinesearchInput;
     RSDsolver->Max_Iteration = 10;
     RSDsolver->Stop_Criterion = GRAD_F;
    RSDsolver->Tolerance = 1e-10;
     RSDsolver->CheckParams();
     RSDsolver->Run();
     delete RSDsolver;
    
    //******************************************OUTPUT*******************************************
    std::ofstream outfile("D:\StudyFiles\Research\Codes\c++ research\GROPT_C\GROPT_C\GROPT_C\TestData.txt");
    for (integer j = 0; j < dim; j++) {
        for (integer i = 0; i < numP; i++) {
            outfile << q1[i+j*numP] << " ";
        }
        outfile << std::endl << std::endl;
    }
    outfile << std::endl;
    for (integer j = 0; j < dim; j++) {
        for (integer i = 0; i < numP; i++) {
            outfile << q2[i+j*numP] << " ";
        }
        outfile << std::endl << std::endl;
    }

    
    printf("testtest\n");
    
    delete [] q1;
    delete [] q2;
}

#endif
