
#include "test/TestElasticCurvesRO.h"

using namespace ROPTLIB;

#if !defined(MATLAB_MEX_FILE) && defined(TESTELASTICCURVESRO)

int main(void)
{
	std::ifstream fdata;
	fdata.open("Curves1020");
	integer NN = 100, num = 1020, n = NN + 1, d = 2;
	double *data = new double[d * NN * num + 2 * d * n];
	double *C1 = data + d * NN * num;
	double *C2 = C1 + d * n;

	for (integer i = 0; i < NN; i++)
	{
		for (integer j = 0; j < num; j++)
		{
			for (integer k = 0; k < 2; k++)
			{
				fdata >> data[k * NN + i + j * 2 * NN];
			}
		}
	}
	// choose a random seed
	unsigned tt = (unsigned)time(NULL);
	tt = 1425718285;
	printf("rand seed:%ud\n", tt);//--
	genrandseed(tt);
	integer idx1 = static_cast<integer> (floor(genrandreal() * 1020)), idx2 = static_cast<integer> (floor(genrandreal() * 1020));
	//idx1 = 0;
	//idx2 = 2;
	printf("idx1:%d, idx2:%d\n", idx1, idx2);
	for (integer i = 0; i < n; i++)
	{
		for (integer j = 0; j < d; j++)
		{
			C1[i + j * n] = data[i + j * NN + idx1 * d * NN];
			C2[i + j * n] = data[i + j * NN + idx2 * d * NN];
		}
	}

	for (integer j = 0; j < d; j++)
	{
		C1[NN + j * n] = C1[0 + j * n];
		C2[NN + j * n] = C2[0 + j * n];
	}

	//double O[] = { cos(0.5), -sin(0.5), sin(0.5), cos(0.5) };
	//char *transn = const_cast<char *> ("n");
	//double one = 1, zero = 0;
	//dgemm_(transn, transn, &n, &d, &d, &one, C1, &n, O, &d, &zero, C2, &n);
	//ForDebug::Print("O first:", O, d, d);//---

	//ForDebug::Print("C1:", C1, n, d);//--
	//ForDebug::Print("C2:", C2, n, d);//--

	//for (integer i = 0; i < n; i++)///---------
	//{
	//	C1[i + 0 * n] = i / (n - 1);
	//	C1[i + 1 * n] = sin(static_cast<double> (i) / (n - 1) / 2);
	//	C2[i + 0 * n] = i / (n - 1);
	//	C2[i + 1 * n] = sin(static_cast<double> (i) / (n - 1) / 2 + PI / 6);
	//}//-------

	integer numofmanis = 3;
	integer numofmani1 = 1;
	integer numofmani2 = 1;
	integer numofmani3 = 1;
	L2SphereVariable *FNSV = new L2SphereVariable(n);
	OrthGroupVariable *OGV = new OrthGroupVariable(d);
	EucVariable *EucV = new EucVariable(1);
	ProductElement *Xopt = new ProductElement(numofmanis, FNSV, numofmani1, OGV, numofmani2, EucV, numofmani3);

	bool rotated = true;
	bool isclosed = true;
	bool onlyDP = false;
	bool swap;
	integer autoselectC = 1;
	double *fopts = new double[10];
	double *comtime = fopts + 5;
	integer ns, lms;

	DriverElasticCurvesRO(C1, C2, d, n, 0.01, rotated, isclosed, onlyDP, 4, "LRBFGS", autoselectC, Xopt, swap, fopts, comtime, ns, lms);

	delete[] fopts;
	delete[] data;
	delete FNSV;
	delete OGV;
	delete EucV;
	delete Xopt;

#ifdef _WIN64
#ifdef _DEBUG
	_CrtDumpMemoryLeaks();
#endif
#endif
	return 0;
};

#endif

#ifdef MATLAB_MEX_FILE

//#define TESTELASTICCURVESRO

std::map<integer *, integer> *CheckMemoryDeleted;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 9)
	{
		mexErrMsgTxt("The number of arguments should be at least nine.\n");
	}
	double *C1, *C2;
	double w = 0;
	C1 = mxGetPr(prhs[0]);
	C2 = mxGetPr(prhs[1]);
	/* dimensions of input matrices */
	integer d, n, rotated, isclosed, onlyDP, skipm, autoselectC;
	n = mxGetM(prhs[0]);
	d = mxGetN(prhs[0]);

	printf("(n, d):%d,%d\n", n, d);

	if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != d)
	{
		mexErrMsgTxt("The size of matrix C2 does not match the size of C1.\n");
	}
	w = mxGetScalar(prhs[2]);
	rotated = static_cast<integer> (mxGetScalar(prhs[3]));
	isclosed = static_cast<integer> (mxGetScalar(prhs[4]));
	onlyDP = static_cast<integer> (mxGetScalar(prhs[5]));
	skipm = static_cast<integer> (mxGetScalar(prhs[6]));
	char methodname[30] = "";
	mxGetString(prhs[7], methodname, 30);
	autoselectC = static_cast<integer> (mxGetScalar(prhs[8]));
	double *Q1 = nullptr, *Q2 = nullptr;
	if(nrhs == 11)
	{
		Q1 = mxGetPr(prhs[9]);
		Q2 = mxGetPr(prhs[10]);
		if (mxGetM(prhs[9]) != n || mxGetN(prhs[9]) != d)
		{
			mexErrMsgTxt("The size of matrix Q1 does not match the size of C1.\n");
		}
		if (mxGetM(prhs[10]) != n || mxGetN(prhs[10]) != d)
		{
			mexErrMsgTxt("The size of matrix Q2 does not match the size of C1.\n");
		}
	}

	genrandseed(0);

	CheckMemoryDeleted = new std::map<integer *, integer>;

	integer numofmanis = 3;
	integer numofmani1 = 1;
	integer numofmani2 = 1;
	integer numofmani3 = 1;
	L2SphereVariable FNSV(n);
	OrthGroupVariable OGV(d);
	EucVariable EucV(1);
	ProductElement *Xopt = new ProductElement(numofmanis, &FNSV, numofmani1, &OGV, numofmani2, &EucV, numofmani3);

    bool swap;
	plhs[2] = mxCreateDoubleMatrix(5, 1, mxREAL);
	plhs[3] = mxCreateDoubleMatrix(5, 1, mxREAL);
	double *fopts = mxGetPr(plhs[2]), *comtime = mxGetPr(plhs[3]);
	integer ns, lms;

	DriverElasticCurvesRO(C1, C2, d, n, w, rotated != 0, isclosed != 0, onlyDP != 0, skipm, methodname,
		autoselectC, Xopt, swap, fopts, comtime, ns, lms, Q1, Q2);

	/*create output matrix*/
	integer sizex = n + d * d + 1;
	plhs[0] = mxCreateDoubleMatrix(sizex, 1, mxREAL);
	double *opt = mxGetPr(plhs[0]);
	plhs[1] = mxCreateDoubleScalar(static_cast<double> (swap));
	plhs[4] = mxCreateDoubleScalar(static_cast<double> (ns));
	plhs[5] = mxCreateDoubleScalar(static_cast<double> (lms));

	const double *Xoptptr = Xopt->ObtainReadData();
	integer inc = 1;
	dcopy_(&sizex, const_cast<double *> (Xoptptr), &inc, opt, &inc);

	delete Xopt;
	
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			printf("Global address: %p, sharedtimes: %d\n", iter->first, iter->second);
	}
	delete CheckMemoryDeleted;
	return;
}

#endif
