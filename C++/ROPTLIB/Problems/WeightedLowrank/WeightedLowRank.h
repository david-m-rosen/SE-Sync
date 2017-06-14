/*
This file defines the class for the problem
min_{X \in R_r^{m by n}} \|X - A\|_W^2, 
where R_r{m by n} is the set of m by n matrices with rank r, A is a given m by n matrix, W is a mn by mn symmetric
positive definite matrix and \|M\|_W^2 is the W weighted norm, i.e., \|M\|_W^2 = vec(M)^T W vec(M) and vec(M) is the vector given by
vectorizing the matrix M.
This problem is used in [ZHGVA2015]
	[ZHGVA2015]: Guifang Zhou, Wen Huang, Kyle A. Gallivan, Paul Van Dooren, Pierre-Antoine Absil. Rank-Constrained Optimization: A Riemannian Manifold Approach,
		In Proceeding of European Symposium on Artificial Neural Networks, Computational Intelligence and Machine Learning, 2015.

Problem --> WeightedLowRank

---- WH
*/

#ifndef WEIGHTEDLOWRANK_H
#define WEIGHTEDLOWRANK_H

#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Manifolds/Element.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"
#include "Manifolds/LowRank/LowRank.h"

/*Define the namespace*/
namespace ROPTLIB{

	class WeightedLowRank : public Problem{
	public:
		WeightedLowRank(double *inA, double *inW, integer inm, integer inn, integer inr);
		virtual ~WeightedLowRank();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		//virtual void RieGrad(Variable *x, Vector *gf) const;

		///*This function is not done yet*/
		//virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		double *A;
		double *W;
		integer m;
		integer n;
		integer r;

	};
}; /*end of ROPTLIB namespace*/
#endif
