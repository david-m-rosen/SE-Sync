/*
This file defines the class for the problem \min_{x \in R^dim} \sum_{i = 1}^num W(i) \|x - D(:, i)\|_2^2

Problem --> EucFrechetMean

---- WH
*/

#ifndef EUCFRECHETMEAN_H
#define EUCFRECHETMEAN_H

#include "Manifolds/Euclidean/Euclidean.h"
#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Problems/Problem.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucFrechetMean : public Problem{
	public:
		EucFrechetMean(double *W, double *D, integer num, integer dim);
		virtual ~EucFrechetMean();
		virtual double f(Variable *x) const;
		virtual void Grad(Variable *x, Vector *gf) const;
		virtual void HessianEta(Variable *x, Vector *etax, Vector *xix) const;

		double *Weights; /* length num */
		double *Data;    /* length num * dim */
		integer Num;
		integer Dim;
	};
}; /*end of ROPTLIB namespace*/
#endif
