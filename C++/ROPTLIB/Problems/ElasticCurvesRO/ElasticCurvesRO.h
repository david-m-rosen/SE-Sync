/*
This file defines the class for the problem in [WGSA2015]
	[WGSA2015]: Wen Huang and K. A. Gallivan and Anuj Srivastava and P.-A. Absil, Riemannian Optimization for Registration of Curves in Elastic Shape Analysis
		Journal of Mathematical Imaging and Vision, doi = "10.1007/s10851-015-0606-8", 2015
Find best reparameterization (rotation) for functions, open and closed curves in R^d
\int_D \|Oq1(t) - q2(\int_0^t l^4(s) ds + m mod 1) l^2(t)\|_2^2 dt + w B(l)
the domain D = [0, 1] for open curves, D is a circle for closed curves.

Problem --> ElasticCurvesRO

---- WH
*/

#ifndef ELASTICCURVESRO_H
#define ELASTICCURVESRO_H

//#define DEBUG_CLIENTBLOCK   new( _CLIENT_BLOCK, __FILE__, __LINE__)

#include "Manifolds/L2Sphere/L2Sphere.h"
#include "Manifolds/L2Sphere/L2SphereVariable.h"
#include "Manifolds/L2Sphere/L2SphereVector.h"
#include "Manifolds/Euclidean/Euclidean.h"
#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"
#include "Manifolds/OrthGroup/OrthGroup.h"
#include "Manifolds/OrthGroup/OrthGroupVariable.h"
#include "Manifolds/OrthGroup/OrthGroupVector.h"
#include "Manifolds/ProductElement.h"
#include "Manifolds/ProductManifold.h"

#include "Others/def.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/Spline.h"

/*Define the namespace*/
namespace ROPTLIB{

	class ElasticCurvesRO : public Problem{
	public:
		/* q1 and q2 are the SRVF of two curves. The curves is in R^d and represented by n points. w is the weight of the barrier function.
		rotated denotes whether the rotation is applied. isclosed denotes whether the curve is closed or not.
		If isclosed is true, then make sure the first and last points of C1 and C2 must be the same.*/
		ElasticCurvesRO(double *inq1, double *inq2, integer ind, integer inn, double inw, bool inrotated, bool inisclosed);

		/*Destructor*/
		virtual ~ElasticCurvesRO();

		//virtual void SetDomain(Manifold *inDomain);
		//	Element **GetInitialXandSetDomain(integer &length);

		/*Cost function*/
		virtual double f(Variable *x) const;

		/*Euclidean gradient*/
		virtual void EucGrad(Variable *x, Vector *egf) const;

		/*Action of the Euclidean Hessian, i.e., exix = EucHess f(x)[etax]*/
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		/*==========The following functions are tools for manipulating the SRV functions=============*/

		/*result(i) = q1(i, :) * q2(i, :)^T for all i */
		static void PointwiseInnerProd(const double *q1, const double *q2, integer d, integer n, double *result);

		/*result(i, j) = q1(i, j) * l(i) for all i \in [0, n - 1] and j \in [0, d - 1] */
		static void PointwiseQProdl(const double *q1, const double *l, integer d, integer n, double *result);

		/*result(i) = l1(i) * l2(i), for all i */
		static void PointwiseProd(const double *l1, const double *l2, integer n, double *result);

		/*Composite trapezoidal rule is used to compute the cumulative integral */
		static void CumTrapz(const double *l, integer n, double intv, double *result);

		/*Composite trapezoidal rule is used to compute the integral */
		static double Trapz(const double *l, integer n, double intv);

		mutable double w;			/* coefficient of penlty term*/
	private:

		double *q1;			/* q1 is represented by n by d matrices*/
		double *q2_coefs;	/* q2 is represented by cubic splines, which represented by (n - 1) by 4 d matrices */
		double *dq2_coefs;	/* The first derivative of the cubic splines, i.e., (n - 1) by 3 d matrices */
		double *ddq2_coefs;	/* The second derivative of the cubic splines, i.e., (n - 1) by 2 d matrices */
		mutable integer n;	/* the number of points for representing a function/curve */
		mutable integer d;	/* the dimension of the curves lie */
		bool rotated;		/* involve rotation or not */
		bool isclosed;		/* for closed curves or open curves*/
	};

}; /*end of ROPTLIB namespace*/
#endif // ELASTICCURVESRO_H
