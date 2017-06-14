/*
This file defines the class for the problem 
\min_{diag(X^T X) = I_r} \mu \|X^T B B^T X - D^2\|_F^2
\|X\|_1
X \in R^{p \times r}, B \in R^{p \times n}. X is a loading matrix. p > n > r.
When proximal update in Oblique class is used, we essentially solve
\min_{diag(X^T X) = I_r} ||X||_1 \mu \|X^T B B^T X - D^2\|_F^2
which is a sparse PCA problem.

Problem --> ObliqueSparsePCA

---- WH
*/

#ifndef OBLIQUESPARSEPCA_H
#define OBLIQUESPARSEPCA_H


#include "Manifolds/Oblique/Oblique.h"
#include "Manifolds/Oblique/ObliqueVariable.h"
#include "Manifolds/Oblique/ObliqueVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ObliqueSparsePCA : public Problem{
	public:
		ObliqueSparsePCA(double *inB, double *inDsq, double inmu, integer inp, integer inn, integer inr);
		virtual ~ObliqueSparsePCA();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *B;
		double *Dsq;
		double mu;
		integer n;
		integer p;
		integer r;
	};

}; /*end of ROPTLIB namespace*/

#endif // end of OBLIQUETESTSPARSEPCA_H