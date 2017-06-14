/*
This file defines the class for problem from julia function handles. This is the interface for julia users.

Problem

---- WH
*/

#ifndef JULIAPROBLEM_H
#define JULIAPROBLEM_H

#include "Problems/Problem.h"
#include <cstring>
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

#ifdef DRIVERJULIAPROB

#include "julia.h"

/*Define the namespace*/
namespace ROPTLIB{

    class juliaProblem : public Problem{
	public:
        /*Construct a julia Problem with function handles of cost function, Euclidean gradient and action of
        Euclidean of Hessian from julia*/
        juliaProblem(jl_function_t *inf, jl_function_t *ingf, jl_function_t *inHess);

		/*Destructor*/
        virtual ~juliaProblem();

        /*call the Julia function handle of the cost function*/
		virtual double f(Variable *x) const;

        /*call the Julia function handle of the Euclidean gradient*/
		virtual void EucGrad(Variable *x, Vector *egf) const;

        /*call the Julia function handle of the action of the Euclidean Hessian*/
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	protected:
        jl_function_t *jl_f; /*Julia function handle of the cost function*/
        jl_function_t *jl_gf; /*Julia function handle of the Euclidean gradient*/
        jl_function_t *jl_Hess; /*Julia function handle of the action of the Euclidean Hessian.*/
	};
}; /*end of ROPTLIB namespace*/


#endif // end of DRIVERJULIAPROB

#endif // end of MEXPROBLEM_H
