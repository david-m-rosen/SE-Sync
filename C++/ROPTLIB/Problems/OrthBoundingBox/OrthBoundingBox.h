/*
This file defines the class for the problem
min_X = V(XE) (e_{1,max} - e_{1, min}) (e_{2, max} - e_{2, min}) ... (e_{n, max} - e_{n, min}),
where X \in O_d, E \in R^{d \times n}, e_{i, max} and e_{i, min} represent the max and min entry
of the i-th row of $XE$, respectively.

Problem --> OrthBoundingBox

---- WH
*/

#ifndef ORTHBOUNDINGBOX_H
#define ORTHBOUNDINGBOX_H

#include "Manifolds/OrthGroup/OrthGroup.h"
#include "Manifolds/OrthGroup/OrthGroupVariable.h"
#include "Manifolds/OrthGroup/OrthGroupVector.h"
#include "Problems/Problem.h"
#include "Manifolds/SharedSpace.h"
#include "Others/def.h"
#include "Others/MyMatrix.h"

/*Define the namespace*/
namespace ROPTLIB{

	class OrthBoundingBox : public Problem{
	public:
		OrthBoundingBox(double *inE, integer ind, integer inn);
		virtual ~OrthBoundingBox();
		virtual double f(Variable *x) const;

		virtual void EucGrad(Variable *x, Vector *egf) const;
		virtual void EucHessianEta(Variable *x, Vector *etax, Vector *exix) const;

		void AllDirectionDerivative(const Variable *x, Vector **dgfs, integer UpperBound, double threshold, integer &idxgf);

		void RecursiveDirDeri(const Variable *x, integer *maxidx, integer *minidx, double *maxv,
			double *minv, integer idx, double volumn, bool findmax, integer &idxgf, Vector **dgfs, integer UpperBound, double *XE, double threshold);

		/*TODO, This can be done by checking whether 0 is in a linear combination of values in each dimension or not.*/
		bool IsStationaryPt(const Variable *x){ return false; };

		double *E;
		integer d;
		integer n;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of ORTHBOUNDINGBOX_H
