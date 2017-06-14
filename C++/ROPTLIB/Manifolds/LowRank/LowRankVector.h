/*
This file defines the class of a point on the tangent space of low-rank manifold R_r^{m times n}

SmartSpace --> ProductElement --> LowRankVector

---- WH
*/

#ifndef LOWRANKVECTOR_H
#define LOWRANKVECTOR_H

#include "Manifolds/ProductElement.h"
#include "Manifolds/Grassmann/GrassVector.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Others/SparseBLAS/blas_sparse.h"

/*Define the namespace*/
namespace ROPTLIB{

	class LowRankVector : public ProductElement{
	public:
		/*Construct an empty vector on the R_r^{m \times n} with only size information.
		The representation of the tangent vector is
		R^{Ur times Uc} times R^{Drc times Drc} times R^{Vr times Vc} */
		LowRankVector(integer Ur, integer Uc, integer Drc, integer Vr, integer Vc);

		/*Destruct by deleting variables
		For low rank tangent vector, the sparse matrix need be deleted by calling
		the sparse blas function: BLAS_usds().*/
		virtual ~LowRankVector(void);

		/*Create an object of LowRankVector with same size as this LowRankVector.*/
		virtual LowRankVector *ConstructEmpty(void) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of OBLIQUEVECTOR_H
