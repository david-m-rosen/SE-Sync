/*
This file defines a function converting complex numbers from
BLAS format to Matlab format

---- WH
*/

#ifdef MATLAB_MEX_FILE

#include "Others/def.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs != 1)
	{
		mexErrMsgTxt("The number of arguments should be one.\n");
	}
	double *B = mxGetPr(prhs[0]);
	/* dimensions of input matrices */
	integer p, n, two = 2, one = 1;
	n = mxGetM(prhs[0]);
	p = mxGetN(prhs[0]);
	integer length = n * p / 2;

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n / 2, p, mxCOMPLEX);
	double *R = mxGetPr(plhs[0]);
	double *I = mxGetPi(plhs[0]);

	dcopy_(&length, B, &two, R, &one);
	dcopy_(&length, B + 1, &two, I, &one);

	return;
}

#endif
