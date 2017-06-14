/*
This file defines the class for the problem
min_{X \in R_r^{m by n}} 0.5 \|P_{\Omega}(X) - P_{\Omega}(A)\|_F^2, 
where R_r{m by n} is the set of m by n matrices with rank r, \Omega is a index set, A_{ij}, (i, j) \in \Omega
are given.

Problem --> LRMatrixCompletion

---- WH
*/

#ifndef LRMATRIXCOMPLETION_H
#define LRMATRIXCOMPLETION_H

#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/SparseBLAS/blas_sparse.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"
#include "Manifolds/LowRank/LowRank.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LRMatrixCompletion : public Problem{
	public:
		LRMatrixCompletion(integer *inir, integer *injc, double *inV, integer innz, integer inm, integer inn, integer inr);

		virtual ~LRMatrixCompletion();

		/*0.5 \|P_omaga(X) - P_omega(A)\|_F^2*/
		virtual double f(Variable *x) const;

		/*P_omaga(X) - P_omega(A)*/
		virtual void EucGrad(Variable *x, Vector *egf) const;

		/*P_omaga(etax)*/
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		/*compute result = P_Omega(U D V^T), where the indices of Omega are stored in inir and injc.*/
		static void ProjecOmegaUDVT(const double *U, const double *D, const double *V, integer inm, integer inn, integer inr, integer *inir, integer *injc, integer nz, double *result);

		integer *ir;
		integer *jc;
		double *V;
		integer nz;
		integer m;
		integer n;
		integer r;

	};
}; /*end of ROPTLIB namespace*/
#endif
