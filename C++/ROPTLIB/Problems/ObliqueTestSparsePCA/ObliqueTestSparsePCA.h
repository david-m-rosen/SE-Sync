/*
This file defines the class for the problem 
\min_{diag(X^T X) = I_r} \|X\|_1 + \mu \|X^T B B^T X - D^2\|_F^2
\|X\|_1 is approximated by a smooth function. (|x| \aprox \sqrt(x^2 + eps^2) - eps)
X \in R^{p \times r}, B \in R^{p \times n}. X is a loading matrix. p > n > r

Problem --> ObliqueTestSparsePCA

---- WH
*/

#ifndef OBLIQUETESTSPARSEPCA_H
#define OBLIQUETESTSPARSEPCA_H


#include "Manifolds/Oblique/Oblique.h"
#include "Manifolds/Oblique/ObliqueVariable.h"
#include "Manifolds/Oblique/ObliqueVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ObliqueTestSparsePCA : public Problem{
	public:
		ObliqueTestSparsePCA(double *inB, double *inDsq, double inmu, double ineps, integer inp, integer inn, integer inr);
		virtual ~ObliqueTestSparsePCA();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *B;
		double *Dsq;
		double mu;
		double epsilon;
		integer n;
		integer p;
		integer r;
	};

}; /*end of ROPTLIB namespace*/

#endif // end of OBLIQUETESTSPARSEPCA_H