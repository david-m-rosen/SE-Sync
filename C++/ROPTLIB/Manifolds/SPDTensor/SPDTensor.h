/*
This file defines the class for the product of symmetric positive definite matrices manifold (SPDTensor).
The affine invariant metric is used. Intrinsic representation of tangent vectors are used.
See details in: 
	[CS15] Anoop Cherian and Suvrit Sra. "Riemannian Dictionary Learning and Sparse Coding for Positive
	Definite Matrices"

Manifold --> SPDTensor

---- WH
*/
#ifndef SPDTENSOR_H
#define SPDTENSOR_H

#include "Manifolds/SPDTensor/SPDTVariable.h"
#include "Manifolds/SPDTensor/SPDTVector.h"
#include "Manifolds/SPDManifold/SPDManifold.h"
#include "Manifolds/ProductManifold.h"
#include "Others/MyMatrix.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SPDTensor : public ProductManifold{
	public:
		/*Construct the SPDTensor manifold*/
		SPDTensor(integer indim, integer innum);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~SPDTensor(void);

		/*X_i = L_i L_i^T, the L is attached to X as a temporary data. The subscript i denotes the i-th slice*/
		void CholeskyRepresentation(Variable *x) const;
	protected:
		integer dim;
		integer num;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of SPDTENSOR_H
