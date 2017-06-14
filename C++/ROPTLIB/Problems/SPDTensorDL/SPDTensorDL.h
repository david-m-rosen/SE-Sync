/*
This file defines the class for the problem: 
the Dictionary Learning problem on the manifold defined in SPDTensor.h and SPDTensor.cpp.
min_{\mathbf{X}} 1/2 * \sum_{j = 1}^N \|log (A_j^{-1/2} (\mathbf{X} alpha_j) A_j^{-1/2})\|_F^2 + lambda_B \sum_{i = 1}^num Tr(X_i),
where X_i is the i-th slice of the tensor \mathbf{X} \in R^{dim \times dim \times num}, alpha_j \in R^{num} is positive
and sparse.
See details in Section IV.A in [CS15].
	[CS15] Anoop Cherian and Suvrit Sra. "Riemannian Dictionary Learning and Sparse Coding for Positive
	Definite Matrices"

Problem --> SPDTensorDL

---- WH
*/

#ifndef SPDTENSORDL_H
#define SPDTENSORDL_H

#include "Manifolds/SPDTensor/SPDTensor.h"
#include "Manifolds/SPDTensor/SPDTVariable.h"
#include "Manifolds/SPDTensor/SPDTVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDTensorDL : public Problem{
	public:
		/*The SPD matrices atoms are stored as their Cholesky matrix, i.e., Li Li^T.
		inLs is a dim by dim by num array of double numbers. alpha is an num by N array
		of double numbers.*/
		SPDTensorDL(double *inLs, double *inalpha, integer indim, integer inN, integer innum, double inlambdaX);

		virtual ~SPDTensorDL();
		
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *gf) const;

		///*Riemannian action of the Hessian has not been done yet*/
		virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		double *Ls;
		double *alpha;
		double lambdaX;
		integer dim;
		integer N;
		integer num;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of STIEBROCKETT_H
