/*
This file defines the class for the problem: 
the Sparse Coding problem on the first quadrant defined in EucPositive.h and EucPositive.cpp.
min_{x >= 0} 1/2 * \|log (A^{-1/2} (\mathbf{B} a) A^{-1/2})\|_F^2 + lambda_a \sum_{i = 1}^num x_i,
where x = (x_1, \ldots, x_num)^T, \mathbf{B} \in R^{dim \times dim \times num} is a tensor,
\mathbf{B} x = \sum B_i a_i and B_i is i-th slice of \mathbf{B}.
See details in Section IV.B in [CS15].
	[CS15] Anoop Cherian and Suvrit Sra. "Riemannian Dictionary Learning and Sparse Coding for Positive
	Definite Matrices"

Problem --> EucPosSpCd

---- WH
*/

#ifndef EUCPOSSPCD_H
#define EUCPOSSPCD_H

#include "Manifolds/EucPositive/EucPositive.h"
#include "Manifolds/EucPositive/EucPosVariable.h"
#include "Manifolds/EucPositive/EucPosVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class EucPosSpCd : public Problem{
	public:
		/*The SPD data points are stored as their Cholesky matrices, i.e., L L^T.
		inLs is a dim by dim by N array of double numbers. inB is an dim by dim by num array
		of double numbers.*/
		EucPosSpCd(double *inLs, double *inB, double inlambdaa, integer indim, integer innum, integer inN);

		virtual ~EucPosSpCd();
		
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *gf) const;

		///*Riemannian action of the Hessian has not been done yet*/
		virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		double *Ls;
		double *B;
		double lambdaa;
		integer dim;
		integer num;
		integer N;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of STIEBROCKETT_H
