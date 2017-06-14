/*
This file defines a function converting complex numbers from
Matlab format to BLAS format

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
	double *R = mxGetPr(prhs[0]);
	double *I = mxGetPi(prhs[0]);
	/* dimensions of input matrices */
	integer p, n, two = 2, one = 1;
	n = mxGetM(prhs[0]);
	p = mxGetN(prhs[0]);
	integer length = n * p;

	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n * 2, p, mxREAL);
	double *ptr = mxGetPr(plhs[0]);

	dcopy_(&length, R, &one, ptr, &two);
	dcopy_(&length, I, &one, ptr + 1, &two);

	return;
}

#endif
