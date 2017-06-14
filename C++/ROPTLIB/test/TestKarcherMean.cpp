//
//  TestKarcherMean.cpp
//  Updated
//
//  Created by Yaqing You on 5/7/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#include "test/TestKarcherMean.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTKARCHERMEAN)

int main(void)
{
    testKarcherMean();
    return 0;
}

void testKarcherMean()
{
    integer numP = 101, dim = 2, numC = 11, numS = 1;
    double *Qs = new double[numP*dim*numS];
    integer NXD = numP*dim;
    std::ifstream infile("/Users/YaqingYou/Documents/Research/Data_Qs.txt");
    if (! infile){
        printf("File did not open\n");
    }
    else{
        for (integer k = 0; k < numS; k++) {
            for (integer i = 0; i < numP; i++) {
                for (integer j = 0; j < dim; j++) {
                    infile >> Qs[i+j*numP+k*numP*dim];
                }
            }
        }
    }

	ForDebug::Print("Qs:", Qs, numP, dim);//---
 //   
    
        //ShapeVariable *shape = new ShapeVariable(numP, dim);
        ElasticShape *Domain = new ElasticShape(numP, dim);
        ShapeVariable *Xinitial = new ShapeVariable(numP, dim);
        double *initial = Xinitial->ObtainWriteEntireData();
    
    std::ifstream infile_2("/Users/YaqingYou/Documents/Research/Data_initial.txt");
    if (! infile_2) {
        printf("File did not open\n");
    }
    else{
        for (integer i = 0; i < numP; i++) {
            for (integer j = 0; j < dim; j++) {
                infile_2 >> initial[i+j*numP];
            }
        }
	}
    
    ForDebug::Print("initial:", initial, numP, dim);//---

    //ForDebug::Print("Qs", Qs, numP, dim);
    //ForDebug::Print("Xinitial", initial, numP, dim);
    KarcherMean Prob(Qs, numP, dim, numC, numS);
    Prob.SetDomain(Domain);
    
    RSD *RSDsolver = new RSD(&Prob, Xinitial);
    RSDsolver->Debug = ITERRESULT;
    RSDsolver->Max_Iteration = 50;
    //RSDsolver->CheckParams();
    RSDsolver->Run();
    const double *opt = RSDsolver->GetXopt()->ObtainReadData();
    std::ofstream outfile("/Users/YaqingYou/Documents/Research/Data_results.txt");
    for (integer j =  0; j < numP; j++) {
        for (integer d = 0; d < dim; d++) {
            outfile << opt[j+d*numP] << ' ';
        }
        outfile << '\n';
    }
    delete RSDsolver;
    
    delete Domain;
    delete Xinitial;
    
    delete [] Qs;

}

#endif

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 6)
	{
		mexErrMsgTxt("The number of arguments should be .\n");
	}
	double *Qs, *initial_in;
	integer numP, dim, numS, numC;

	Qs = mxGetPr(prhs[0]);
	initial_in = mxGetPr(prhs[1]);
	numP = static_cast<integer>(mxGetScalar(prhs[2]));
	dim = static_cast<integer>(mxGetScalar(prhs[3]));
	numS = static_cast<integer>(mxGetScalar(prhs[4]));
	numC = static_cast<integer>(mxGetScalar(prhs[5]));

	integer NXD = numP*dim;

	printf("numP %d dim %d numS %d numC %d", numP, dim, numS, numC);
	printf("\n");
// 	printf("Qs:=====================================\n");
// 	printf("%f ", Qs[0]);
// 	printf("%f ", Qs[1]);
// 	printf("%f ", Qs[2]);

// 	printf("\n");
// 	printf("qinitial_in:=====================================\n");
// 	for (integer i = 0; i < NXD; i++)
// 	{
// 		printf("%f ", initial_in[i]);
// 	}



	ElasticShape *Domain = new ElasticShape(numP, dim);
	ShapeVariable *Xinitial = new ShapeVariable(numP, dim);
	double *initial = Xinitial->ObtainWriteEntireData();
	dcopy_(&NXD, initial_in, &GLOBAL::IONE, initial, &GLOBAL::IONE);

	printf("initial:=====================================\n");
	for (integer i = 0; i < NXD; i++)
	{
		printf("%f ", initial[i]);
	}


	KarcherMean Prob(Qs, numP, dim, numC, numS);

	Prob.SetDomain(Domain);
	RSD *RSDsolver = new RSD(&Prob, Xinitial);
	RSDsolver->Debug = ITERRESULT;
	RSDsolver->Max_Iteration = 20;
    RSDsolver->Stop_Criterion = FUN_REL;
    RSDsolver->Tolerance = 5e-5;
	//RSDsolver->CheckParams();
	RSDsolver->Run();

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(numP, dim, mxREAL);
	double *opt = mxGetPr(plhs[0]);
	const double *oopt = RSDsolver->GetXopt()->ObtainReadData();
	dcopy_(&NXD, const_cast<double *>(oopt), &GLOBAL::IONE, opt, &GLOBAL::IONE);


	delete RSDsolver;
	delete Domain;
	delete Xinitial;

	return;
}

#endif
