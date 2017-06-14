/*
This file defines the class for the oblique manifold Ob(n, num) = \{X \in R^{n times num} | diag(X^T X) = I_{num} \}.
Note that we don't force X to have full-column rank. It follows that X is actually just a product of unit spheres.

Manifold --> ProductManifold --> Oblique

---- WH
*/

#ifndef OBLIQUE_H
#define OBLIQUE_H

#include "Manifolds/ProductManifold.h"
#include "Manifolds/Sphere/Sphere.h"
#include "Manifolds/Oblique/ObliqueVariable.h"
#include "Manifolds/Oblique/ObliqueVector.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Oblique : public ProductManifold{
	public:
		/*Construct an Oblique manifold, which is a product of spheres.
		the number of sphere is num and the unit sphere is in R^n.*/
		Oblique(integer n, integer num);

		/*Delete each component manifold*/
		~Oblique();

		/* Default one: choose qf, parallelization and intrinsic approach and use householder reflections*/
		virtual void ChooseObliqueParamsSet1();

		/* choose exponential map, parallel translation and extrinsic approach and no householder reflections
		Even though the Householder reflections are not used, the locking condition is satisfied.*/
		void ChooseObliqueParamsSet2();

		/* choose qf, parallel translation and extrinsic approach and no householder reflections
		The locking conidition is not satisfied*/
		void ChooseObliqueParamsSet3();

		/* choose qf, parallel translation and extrinsic approach and no householder reflections
		Beta \neq 1 is used and the locking conidition is satisfied*/
		void ChooseObliqueParamsSet4();

		/* choose proximal mapping, parallelization and intrinsic approach and no householder reflections
		the locking conidition is not satisfied. This is used to minimize f(X) + ||X||_1, where f(X) is
		the original cost function defined on a manifold*/
		void ChooseObliqueParamsSet5();

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);
	};
}; /*end of ROPTLIB namespace*/
#endif // end of OBLIQUE_H
