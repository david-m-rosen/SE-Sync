/*
This file defines the class for the problem (details in (12.1.1) in Wen Huang's Thesis)
\min_{X \in \St(p, n)} - \Sum_{i=1}^N \| diag(X^T C_i X^T) \|_F^2,
where St(p, n) is the Stiefel manifold, C_i are given symmetric matricies.

Problem --> StieSoftICA

---- WH
*/

#ifndef STIESOFTICA_H
#define STIESOFTICA_H

#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class StieSoftICA : public Problem{
	public:
		StieSoftICA(double *inCs, integer inn, integer inp, integer inN);
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		double *Cs;
		mutable integer n;
		mutable integer p;
		mutable integer N;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of STIESOFTICA_H
