
#include "Problems/mexProblem.h"

#ifdef MATLAB_MEX_FILE

/*Define the namespace*/
namespace ROPTLIB{

	mexProblem::mexProblem(const mxArray *inf, const mxArray *ingf, const mxArray *inHess)
	{
		mxf = inf;
		mxgf = ingf;
		mxHess = inHess;
	};

	mexProblem::~mexProblem()
	{
	};

	double mexProblem::f(Variable *x) const
	{
		mxArray *Xmx;
		ObtainMxArrayFromElement(Xmx, x);
		mxArray *lhs[2], *rhs[2];
		rhs[0] = const_cast<mxArray *> (mxf);
		rhs[1] = const_cast<mxArray *> (Xmx);
		mexCallMATLAB(2, lhs, 2, rhs, "feval");
		ObtainElementFromMxArray(x, lhs[1]);
		double result = mxGetScalar(lhs[0]);
		mxDestroyArray(Xmx);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
		return result;
	};

	void mexProblem::EucGrad(Variable *x, Vector *egf) const
	{
		mxArray *Xmx;
		ObtainMxArrayFromElement(Xmx, x);
		mxArray *lhs[2], *rhs[2];
		rhs[0] = const_cast<mxArray *> (mxgf);
		rhs[1] = const_cast<mxArray *> (Xmx);
		mexCallMATLAB(2, lhs, 2, rhs, "feval");
		ObtainElementFromMxArray(x, lhs[1]);
		ObtainElementFromMxArray(egf, lhs[0]);
		mxDestroyArray(Xmx);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
	};

	void mexProblem::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		mxArray *Xmx, *Etaxmx;
		ObtainMxArrayFromElement(Xmx, x);
		ObtainMxArrayFromElement(Etaxmx, etax);
		mxArray *lhs[2], *rhs[3];
		rhs[0] = const_cast<mxArray *> (mxHess);
		rhs[1] = const_cast<mxArray *> (Xmx);
		rhs[2] = const_cast<mxArray *> (Etaxmx);
		mexCallMATLAB(2, lhs, 3, rhs, "feval");
		ObtainElementFromMxArray(x, lhs[1]);
		ObtainElementFromMxArray(exix, lhs[0]);
		mxDestroyArray(Xmx);
		mxDestroyArray(Etaxmx);
		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
	};

	void mexProblem::ObtainMxArrayFromElement(mxArray *&Xmx, const Element *X)
	{
		integer sizeoftempdata = X->GetSizeofTempData();
		integer nfields = sizeoftempdata + 1;
		std::string *fnames = new std::string[nfields];
		X->ObtainTempNames(fnames);
		fnames[sizeoftempdata].assign("main");
		char **names = new char *[nfields];
		for (integer i = 0; i < nfields; i++)
		{
			names[i] = new char[30];
			if (fnames[i].size() >= 30)
			{
				mexErrMsgTxt("The lengths of field names should be less than 30!");
			}
			strcpy(names[i], fnames[i].c_str());
		}

		mxArray *tmp;
		const char *name;
		const SharedSpace *Sharedtmp;
		integer lengthtmp, inc = 1;
		const double *Sharedtmpptr;
		double *tmpptr;
		Xmx = mxCreateStructMatrix(1, 1, nfields, const_cast<const char **> (names));
		for (integer i = 0; i < nfields; i++)
		{
			name = mxGetFieldNameByNumber(Xmx, i);
			if (strcmp(name, "main") != 0)
			{
				Sharedtmp = X->ObtainReadTempData(name);
				Sharedtmpptr = Sharedtmp->ObtainReadData();
				lengthtmp = Sharedtmp->Getlength();
			}
			else
			{
				Sharedtmpptr = X->ObtainReadData();
				lengthtmp = X->Getlength();
			}
			tmp = mxCreateDoubleMatrix(lengthtmp, 1, mxREAL);
			tmpptr = mxGetPr(tmp);
			// tmpptr <- Sharedtmpptr, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
			dcopy_(&lengthtmp, const_cast<double *> (Sharedtmpptr), &inc, tmpptr, &inc);
			mxSetFieldByNumber(Xmx, 0, i, tmp);
		}

		for (integer i = 0; i < nfields; i++)
			delete[] names[i];
		delete[] names;
		delete[] fnames;
	};

	void mexProblem::ObtainElementFromMxArray(Element *X, const mxArray *Xmx)
	{
		// copy the main data from mxArray to X
		double *Xmxptr = mxGetPr(GetFieldbyName(Xmx, 0, "main"));
		integer lengthX = X->Getlength();
		integer inc = 1;
		double *Xptr = X->ObtainWriteEntireData();
		dcopy_(&lengthX, Xmxptr, &inc, Xptr, &inc);

		// copy temp data from mxArray to X
		integer nfields = mxGetNumberOfFields(Xmx);
		const char **fnames;
		mxArray *tmp;

		fnames = static_cast<const char **> (mxCalloc(nfields, sizeof(*fnames)));
		for (integer i = 0; i < nfields; i++)
		{
			fnames[i] = mxGetFieldNameByNumber(Xmx, i);
			if (strcmp(fnames[i], "main") != 0)
			{
				tmp = GetFieldbyName(Xmx, 0, fnames[i]);
				double *tmpptr = mxGetPr(tmp);
				integer length = mxGetM(tmp);
				SharedSpace *Sharedtmp = new SharedSpace(1, length);
				double *Sharedtmpptr = Sharedtmp->ObtainWriteEntireData();
				// Sharedtmpptr <- tmpptr, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
				dcopy_(&length, tmpptr, &inc, Sharedtmpptr, &inc);
				X->AddToTempData(fnames[i], Sharedtmp);
			}
		}
	};

	mxArray *mexProblem::GetFieldbyName(const mxArray *S, integer idxstruct, const char *name)
	{
		integer nfields = mxGetNumberOfFields(S);
		for (integer i = 0; i < nfields; i++)
		{
			if (strcmp(mxGetFieldNameByNumber(S, i), name) == 0)
			{
				return mxGetFieldByNumber(S, idxstruct, i);
			}
		}
		return nullptr;
	};
}; /*end of ROPTLIB namespace*/

#endif
