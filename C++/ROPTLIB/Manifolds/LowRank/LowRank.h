/*
This file defines the class for the low-rank manifold R_r^{m times n}, which is represented by Gr(r, m) times R^{r times r} times Gr(r, n).
A tangent vector in T_x R_r^{m times n} can representated by 
etax = dot{U} D V^T + U dot{D} V^T + U D dot{V}^T, where dot{U}^T U = 0, dot{V}^T V = 0.

In the implementation, a tangent vector has 3 representations:
1), extrinsic representation: (dot{U}, dot{D}, dot{V})
2), intrinsic representation: a (mr + nr - r^2) vector
3), Euclidean representation: a m-by-n matrix.
Given a Euclidean gradient, which is usually represented in the third form,
The Euclidean representation is stored in a temporary data with key "EucRep" in
a EMPTYEXTR-type tangent vector.
This class provides a function to convert the Euclidean representation into the extrinsic representation.

The used Riemannian metric is
g(etax, xix) = trace(etax^T xix)
= trace(D^T \dot{U}_1^T \dot{U}_2 D) + \trace(\dot{D}_1^T \dot{D}_2) + \trace(D \dot{V}_1^T \dot{V}_2 D^T),
where etax = dot{U}_1 D V^T + U dot{D}_1 V^T + U D dot{V}_1^T and xix = dot{U}_2 D V^T + U dot{D}_2 V^T + U D dot{V}_2^T.

Manifold --> ProductManifold --> LowRank

---- WH
*/

#ifndef LOWRANK_H
#define LOWRANK_H

#include "Manifolds/ProductManifold.h"
#include "Manifolds/Stiefel/Stiefel.h"
#include "Manifolds/Grassmann/Grassmann.h"
#include "Manifolds/Euclidean/Euclidean.h"
#include "Manifolds/LowRank/LowRankVariable.h"
#include "Manifolds/LowRank/LowRankVector.h"
#include "Others/MyMatrix.h"
#include "Others/SparseBLAS/blas_sparse.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LowRank : public ProductManifold{
	public:
		/*Construct the low rank manifold of m by n matrices with rank r.
		It is represented by St(r, m) times R^{r times r} times St(r, n), i.e.,
		X = U D V^T. U \in St(r, m), D \in R^{r times r} and V \in St(r, n).
		Note that D is not necessary a diagonal matrix.*/
		LowRank(integer m, integer n, integer r);

		/*Delete the manifold by deleting each component.*/
		~LowRank(void);

		/*Riemannian metric*/
		virtual double ExtrMetric(Variable *x, Vector *etax, Vector *xix) const;

		/*Tangent vector etax = \dot(U) D V^T + U \dot{D} V^T + U D \dot{V}^T. Let \dot{U} = U_perp K_U and \dot{V} = V_perp K_V
		It follows that etax = U_perp K_U D V^T + U D K_V^T V_perp^T + U \dot{D} V^T
		The intrinsic representation would be given by vectorizing (A, B, C) := (K_U D, D K_V^T, \dot{D}) .*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Compute the extrinsic approach by given (A, B, C) := (K_U D, D K_V^T, \dot{D}).
		\dot{U} = U_\perp A D^{-1}, \dot{D} = C, \dot{V} = V_\perp B^T D^{-T}.*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*Perform the default retraction of each manifold component*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*Compute the tangent vector result satisfying
		g_y(\mathcal{T}_{R_etax}(xix), xiy) = g_x(xix, result) for all xix \in T_x M,
		where y = R_x(etax), xiy \in T_y M.*/
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Used in "coTangentVectorCan". (TODO: add details in notes?)*/
		virtual void ExtrProjectionStiePerp(Variable *x, Vector *v, Vector *result) const;

		/*Perform the vector transport by differentiated retraction of each individual manifold.*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir) const;
		
		/*etax is in the ambient space R^{m \times r} \times R^{r \times r} \times R^{n \times r}. This function projects etax onto the 
		tangent space of x, i.e., result = P_{T_x M} v, where P is based on the selected Riemannian metric*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*the Riemannian gradient is obtained by projecting the Euclidean onto the tangent space of x.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*The Riemannian action of the Hessian is obtained by Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*Convert the Euclidean representation R^{m \times n} of a tangent vector to the Extrinsic representation.
		The Euclidean representation is attached in the temporary data of "result" with key "EucRep". */
		virtual void EucRepToExtr(Variable *x, Vector *result) const;

		/*Convert the Extrinsic representation of a tangent vector to the Euclidean representation.
		The Euclidean representation is attached in the temporary data of "result" with key "EucRep". */
		virtual void ExtrToEucRep(Variable *x, Vector *result) const;

		/*Construct the sparse matrix using sparse BLAS library*/
		blas_sparse_matrix ConstructSparseMatrix(Vector *result) const;
		
		/*Compute the LU decomposition of x*/
		void LUofDinx(Variable *x) const;
	protected:
		integer m; /*the number of row*/
		integer n; /*the number of column*/
		integer r; /*the rank of the matrix*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of LOWRANK_H
