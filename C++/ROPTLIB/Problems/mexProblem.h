/*
This file defines the class for problem from matlab handles. This is the interface for Matlab users.

Problem

---- WH
*/

#ifndef MEXPROBLEM_H
#define MEXPROBLEM_H

#include "Problems/Problem.h"
#include <cstring>
#include "Others/def.h"

#ifdef MATLAB_MEX_FILE

/*Define the namespace*/
namespace ROPTLIB{

	class mexProblem : public Problem{
	public:
		/*Construct a mex Problem with function handles of cost function, Euclidean gradient and action of
		Euclidean of Hessian from Matlab*/
		mexProblem(const mxArray *inf, const mxArray *ingf, const mxArray *inHess);

		/*Destructor*/
		virtual ~mexProblem();

		/*call the Matlab function handle of the cost function*/
		virtual double f(Variable *x) const;

		/*call the Matlab function handle of the Euclidean gradient*/
		virtual void EucGrad(Variable *x, Vector *egf) const;

		/*call the Matlab function handle of the action of the Euclidean Hessian*/
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		/*This function converts the format of storage in ROPTLIB to the format
		of storage in Matlab*/
		static void ObtainMxArrayFromElement(mxArray *&Xmx, const Element *X);

		/*This function converts the format of storage in Matlab to the format
		of storage in ROPTLIB*/
		static void ObtainElementFromMxArray(Element *X, const mxArray *Xmx);

		/*S is a Matlab structure. This function obtain its field by key = name */
		static mxArray *GetFieldbyName(const mxArray *S, integer idxstruct, const char *name);
	protected:
		const mxArray *mxf; /*Matlab function handle of the cost function*/
		const mxArray *mxgf; /*Matlab function handle of the Euclidean gradient*/
		const mxArray *mxHess; /*Matlab function handle of the action of the Euclidean Hessian.*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of MATLAB_MEX_FILE

#endif // end of MEXPROBLEM_H
