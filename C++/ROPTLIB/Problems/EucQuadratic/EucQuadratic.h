/*
This file defines the class for the problem min_{x \in R^d} x^T A x, where A is a d by d symmetric positive definite matrix

Problem --> EucQuadratic

---- WH
*/

#ifndef EUCQUADRATIC_H
#define EUCQUADRATIC_H

#include "Manifolds/Euclidean/Euclidean.h"
#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucQuadratic : public Problem{
	public:
		EucQuadratic(double *M, integer dim);
		virtual ~EucQuadratic();
		virtual double f(Variable *x) const;
		virtual void Grad(Variable *x, Vector *gf) const;
		virtual void HessianEta(Variable *x, Vector *etax, Vector *xix) const;

		double *A;
		integer Dim;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of EUCQUADRATIC_H
