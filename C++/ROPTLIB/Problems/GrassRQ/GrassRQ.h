/*
This file defines the class for the Rayleigh Quotient problem on the Grassmann manifold
min_X tr(X^T B X), where B is a symmetric matrix, and X \in Gr(p, n).

Problem --> GrassRQ

---- WH
*/

#ifndef GRASSRQ_H
#define GRASSRQ_H

#include "Manifolds/Grassmann/Grassmann.h"
#include "Manifolds/Grassmann/GrassVariable.h"
#include "Manifolds/Grassmann/GrassVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class GrassRQ : public Problem{
	public:
		GrassRQ(double *inB, integer inn, integer inp);
		virtual ~GrassRQ();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *B;
		integer n;
		integer p;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of GRASSRQ_H
