/*
This file defines the class for the problem 
min_{X1 \in St(p, n), X2 \in St(p, n), X3 \in St(q, m)}     tr(X1^T B1 X1 D1) + tr(X2^T B2 X2 D2) + tr(X3^T B3 X3 D3), 
where B1, B2 are n by n symmetric matrices, B3 is m by m symmetric matrix, D1 and D2 are p by p diagonal matrices, and
D3 is q by q diagonal matrix.

Problem --> StieSumBrockett

---- WH
*/

#ifndef STIESUMBROCKETT_H
#define STIESUMBROCKETT_H

#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class StieSumBrockett : public Problem{
	public:
		StieSumBrockett(double *inB1, double *inD1, double *inB2, double *inD2, double *inB3,
			double *inD3, integer inn, integer inp, integer inm, integer inq);
		virtual ~StieSumBrockett();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *B1;
		double *D1;
		double *B2;
		double *D2;
		double *B3;
		double *D3;
		integer n;
		integer p;
		integer m;
		integer q;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of STIESUMBROCKETT_H
