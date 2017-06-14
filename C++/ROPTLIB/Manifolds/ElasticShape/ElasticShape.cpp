//
//  ElasticShape.cpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#include "Manifolds/ElasticShape/ElasticShape.h"

/*Define the namespace*/
namespace ROPTLIB{
    ElasticShape::ElasticShape(integer r, integer c)
    {
        numP = r;
        dim = c;
        IsIntrApproach = false; //tangent vector intrinsic
        HasHHR = false; //householder reflector
        UpdBetaAlone = false;
        name.assign("ElasticShape");
        IntrinsicDim = r * c;
        ExtrinsicDim = r * c;
        EMPTYEXTR = new ShapeVector(r, c);
        EMPTYINTR = new ShapeVector(r, c);
    };

    ElasticShape::~ElasticShape(void)
    {
        delete EMPTYEXTR;
        delete EMPTYINTR;
    }

    double ElasticShape::Metric(Variable *x, Vector *etax, Vector *xix) const
    {
        const double *vec_1 = etax->ObtainReadData();
        const double *vec_2 = xix->ObtainReadData();
        double result;
        result = PreShapeCurves::InnerProd_Q(vec_1, vec_2, numP, dim);
        return result;
    }


    void ElasticShape::Projection(Variable *x, Vector *etax, Vector *result) const
    {
        const double *Curve_x = x->ObtainReadData();
        etax->CopyTo(result);
        //const double *Vec_Etax = etax->ObtainReadData();
        double *ans = result->ObtainWritePartialData();
        /*Project arbitrary points in L2 into tangent space at q*/
        /*void Item_2(const double *q, integer innumP, integer indim, double *w);*/
        PreShapePathStraighten::Item_2(Curve_x, numP, dim, ans);
    }

    void ElasticShape::EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const
    {
        Projection(x, egf, gf);
    }

    void ElasticShape::EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* result, const Problem *prob) const
    {
        exix->CopyTo(result);
    };


    void ElasticShape::Retraction(Variable *x, Vector *etax, Variable *result) const
    {
        const double *Curve_x = x->ObtainReadData();
        const double *direc = etax->ObtainReadData();
        double *Curve_new = result->ObtainWriteEntireData();
        double *Curve_temp = new double[numP*dim];
        for (integer i = 0; i < numP*dim; i++)
        {
            Curve_temp[i] = Curve_x[i] + direc[i];
        }
        PreShapePathStraighten::Item_1(Curve_temp, numP, dim, Curve_new);

        delete[] Curve_temp;
	};

}; /*end of ROPTLIB namespace*/
