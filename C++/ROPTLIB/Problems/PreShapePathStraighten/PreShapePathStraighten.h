/*
 This file defines the class for the problem in [YHGA2015]
 [YHGA2015]: Yaqing You and Wen Huang and Kyle A. Gallivan and P.-A. Absil, A Riemannian Approach for Computing Geodesics in Elastic Shape Analysis, in Proceedings of the 3rd IEEE Global Conference on Signal & Information Processing, 2015
 Find a geodesic between two closed curves in preshape space, using algorithm given in Anuj Srivastava, Eric Klassen, Shantanu H. Joshi and Ian H. Jermyn, Shape Analysis of Elastic Curves in Euclidean Spaces.
 
 Problem -->PreShapePathStraighten
 
 ----- YY
 */


#ifndef PRESHAPEPATHSTRAIGHTEN_H
#define PRESHAPEPATHSTRAIGHTEN_H

#include "Manifolds/PreShapeCurves/PreShapeCurves.h"
#include "Manifolds/PreShapeCurves/PSCVariable.h"
#include "Manifolds/PreShapeCurves/PSCVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Problems/ElasticCurvesRO/ElasticCurvesRO.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"


/*Define the namespace*/
namespace ROPTLIB{
	class PreShapePathStraighten : public Problem{
	public:
		/* numP is the number of points on each curve
		   dim is the dimension of the space where the curve is in
		   numC is the number of curves on each path
		 */
		PreShapePathStraighten(integer innumP, integer indim, integer innumC);
    
		/*Destructor*/
		virtual ~PreShapePathStraighten();
    
		/*cost function*/
		virtual double f(Variable *x) const;

		//virtual void RieGrad(Variable *x, Vector *gf) const;
		//virtual void RieHessianEta(Variable *x, Vector *etax, Vector *xix) const;

		/*Euclidean gradient*/
		virtual void EucGrad(Variable *x, Vector *egf) const;
    
		/*Project a curve onto closed curve space*/
		static void Item_1(const double *q, integer innumP, integer indim, double *q_c);
    
		/*Project arbitrary points in L2 into tangent space at q*/
		static void Item_2(const double *q, integer innumP, integer indim, double *w);
    
		/*Parallel Transport*/
		static void Item_3(const double *w, const double *q1, const double *q2, integer innumP, integer indim, double *wbar);
    
		static double InnerProd_Q(const double *q1, const double *q2, integer innumP, integer indim);
    
		//virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

	private:
		integer numP;
		integer dim;
		integer numC;
	};
}; /*end of ROPTLIB namespace*/

#endif // end of PRESHAPEPATHSTRAIGHTEN_H
