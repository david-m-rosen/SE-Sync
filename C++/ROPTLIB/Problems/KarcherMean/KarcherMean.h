//
//  KarcherMean.hpp
//  Updated
//
//  Created by Yaqing You on 5/4/16.
//  Copyright © 2016 Yaqing You. All rights reserved.
//

#ifndef KarcherMean_h
#define KarcherMean_h

#include "Manifolds/ElasticShape/ElasticShape.h"
#include "Manifolds/ElasticShape/ShapeVariable.h"
#include "Manifolds/ElasticShape/ShapeVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Problems/ShapePathStraighten/ShapePathStraighten.h"

#include <stdio.h>

#include <fstream>

#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "Manifolds/PreShapeCurves/PSCVariable.h"
#include "Manifolds/PreShapeCurves/PSCVector.h"

#include "Manifolds/PreShapeCurves/PreShapeCurves.h"
#include "Problems/PreShapePathStraighten/PreShapePathStraighten.h"

#include "Solvers/RSD.h"
#include "Solvers/RNewton.h"
#include "Solvers/RCG.h"
#include "Solvers/RBroydenFamily.h"
#include "Solvers/RWRBFGS.h"
#include "Solvers/RBFGS.h"
#include "Solvers/LRBFGS.h"

#include "Solvers/SolversTR.h"
#include "Solvers/RTRSD.h"
#include "Solvers/RTRNewton.h"
#include "Solvers/RTRSR1.h"
#include "Solvers/LRTRSR1.h"

#include "Problems/ElasticCurvesRO/DriverElasticCurvesRO.h"

/*Define the namespace*/
namespace ROPTLIB{
    class KarcherMean : public Problem{
    public:
        KarcherMean(double *Curves, integer innumP, integer indim, integer innumC, integer innumS);
        
        virtual ~KarcherMean();
        
        virtual double f(Variable *x) const;
        
        virtual void EucGrad(Variable *x, Vector *egf) const;
        
        void ComputeGeodesic(const double *q1, const double *q2, double &distance, double *w) const;
        
    private:
        integer numC;
        integer numP;
        integer dim;
        integer numS;
        double *Qs;
	};
    double LinesearchInput(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *Prob, const Solvers *solver);
    
}; /*end of ROPTLIB namespace*/

#endif /* KarcherMean_h */
