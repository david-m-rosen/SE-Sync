/*
 This file defines the class for preshape manifold.
 It defines the common properties and features of the manifold
 
 Manifold --> PreShapeCurves
 
 ---- YY
 */

#ifndef PRESHAPECURVES_H
#define PRESHAPECURVES_H

#include "Manifolds/PreShapeCurves/PSCVariable.h"
#include "Manifolds/PreShapeCurves/PSCVector.h"
#include "Manifolds/Manifold.h"
#include "Others/def.h"
#include "Problems/ElasticCurvesRO/ElasticCurvesRO.h"
#include "Problems/PreShapePathStraighten/PreShapePathStraighten.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{
	class PreShapeCurves : public Manifold{
	public:
		/*Construct the preshape manifold and set up parameters*/
		PreShapeCurves(integer r, integer c, integer n);
    
		/*Destructor*/
		virtual ~PreShapeCurves();
    
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* result, const Problem *prob) const;
    
		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;
    
		/*When the metric is Euclidean, the Riemannian gradient is obtained by projecting the Euclidean onto the tangent space of x.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;
    
		static double InnerProd_Q(const double *q1, const double *q2, integer innumP, integer indim);
    
		/*Compute the covariant integral of d\alpha / d\tau, denoted by u*/
		static void CovIntegral(const double *Dalpha, const double *alpha, integer innumC, integer innumP, integer indim, double *u);
    
		/*Compute the backward parallel transport of the vector u(1) along \alpha, denoted by utilde*/
		static void BackTrans(const double *u, const double *alpha, integer innumC, integer innumP, integer indim, double *utilde);
    
		/*Compute the Gradient vector field of the energy E*/
		static void GradVec(const double *utilde, const double *u, integer innumC, integer innumP, integer indim, double *w);

	private:
		integer numP; // number of points
		integer dim;  // dimension
		integer numC; // number of curves on a path
	};
}; /*end of ROPTLIB namespace*/

#endif // end of EUCLIDEAN_H
