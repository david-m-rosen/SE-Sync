/*
This file defines the class for the manifold C_*^{n \times p} / O_p, where C_*^{n \times p} is a n by p full column rank
complex matrix and O_p is a p-by-p unitary group. This manifold is equivalent to the manifold of the set of hermitian positive 
semidefinite matricies with rank fixed (rank p). Details can be found in [HGZ2015].
	[HGZ2015]:W. Huang, Kyle A. Gallivan, and Xiangxiong Zhang. Solving PhaseLift by low rank Riemannian optimization methods for complex semidefinite constraints.
		U.C.Louvain, UCL-INMA-2015.01, 2015.

Manifold --> CpxNStQOrth

---- WH
*/

#ifndef CPXNSTQORTH_H
#define CPXNSTQORTH_H

#include "Manifolds/CpxNStQOrth/CSOVariable.h"
#include "Manifolds/CpxNStQOrth/CSOVector.h"
#include "Manifolds/Manifold.h"
#include "Others/MyMatrix.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class CpxNStQOrth : public Manifold{
	public:
		/*Construct the manifold C_*^{n \times p} / O_p*/
		CpxNStQOrth(integer n, integer p = 1);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~CpxNStQOrth(void);

		/*Euclidean metric*/
		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		/*Call a member function "IntrProjection" or "ExtrProjection" based on member variable "IsIntrApproach"*/
		virtual void Projection(Variable *x, Vector *v, Vector *result) const;

		/*Only support intrinsic representataion for etax. The retraction is result = x + etax.*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*This is not done yet. */
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*The vector transport by differentiated retraction is the same is the vector tranport by projection for this retraction and manifold.*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Only use one*/
		virtual double Beta(Variable *x, Vector *etax) const;

		/*Householder transformations are used to obtain the intrinsic representation of etax*/
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Householder transformations are used to obtain the extrinsic representation of etax*/
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*IntrProjection is identity, i.e., result <-- v*/
		virtual void IntrProjection(Variable *x, Vector *v, Vector *result) const;

		/*ExtrProjection is given in e.g., [HGZ2015, Lemma 7].
			[HGZ2015]:W. Huang, Kyle A. Gallivan, and Xiangxiong Zhang. Solving PhaseLift by low rank Riemannian optimization methods for complex semidefinite constraints.
			U.C.Louvain, UCL-INMA-2015.01, 2015.*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*When the metric is Euclidean, the Riemannian gradient is obtained by projecting the Euclidean
		onto the horizontal space of x. See [HGZ2015, Lemma 9].
		[HGZ2015]:W. Huang, Kyle A. Gallivan, and Xiangxiong Zhang. Solving PhaseLift by low rank Riemannian optimization methods for complex semidefinite constraints.
		U.C.Louvain, UCL-INMA-2015.01, 2015.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*See details in [HGZ2015, Lemma 9].
			[HGZ2015]:W. Huang, Kyle A. Gallivan, and Xiangxiong Zhang. Solving PhaseLift by low rank Riemannian optimization methods for complex semidefinite constraints.
			U.C.Louvain, UCL-INMA-2015.01, 2015.*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;
	protected:
		/*Compute the units vectors in Hourseholder transformations, which can be used in the member functions "ObtainIntr" and "ObtainExtr". */
		void ComputeHHR(Variable *x) const;

		integer n; /*The row*/
		integer p; /*The column*/
	};
}; /*end of ROPTLIB namespace*/

#endif // end of EUCLIDEAN_H
