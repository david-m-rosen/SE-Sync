//
//  ElasticShape.hpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#ifndef ElasticShape_h
#define ElasticShape_h


#include "Manifolds/ElasticShape/ShapeVariable.h"
#include "Manifolds/ElasticShape/ShapeVector.h"
#include "Manifolds/Manifold.h"
#include "Others/def.h"
#include "Problems/PreShapePathStraighten/PreShapePathStraighten.h"
#include "Manifolds/PreShapeCurves/PreShapeCurves.h"

/*Define the namespace*/
namespace ROPTLIB{
	class ElasticShape : public Manifold{
	public:
		/*Construct the elastic shape manifold*/
		ElasticShape(integer r, integer c);

		/*Destructor*/
		virtual ~ElasticShape();

		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		virtual void Projection(Variable *x, Vector *etax, Vector *result) const;

		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* result, const Problem *prob) const;

		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

	private:
		integer numP;
		integer dim;
	};
}; /*end of ROPTLIB namespace*/
#endif /* ElasticShape_h */
