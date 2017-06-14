/*
This file defines the class for the unit sphere S^{n-1} = \{x \in R^{n} | x^T x = 1\}
It defines the common properties and features of the manifold.

Manifold --> Stiefel --> Sphere

---- WH
*/

#ifndef SPHERE_H
#define SPHERE_H

#include "Manifolds/Sphere/SphereVariable.h"
#include "Manifolds/Sphere/SphereVector.h"
#include "Manifolds/Stiefel/Stiefel.h"

/*Define the namespace*/
namespace ROPTLIB{

	class Sphere : public Stiefel{
	public:
		/*Construct the unit sphere S^{n-1}*/
		Sphere(integer n);

		/*Delete the sphere*/
		virtual ~Sphere();

		// choose qf retraction, parallelization and intrinsic approach and no householder reflections
		void ChooseSphereParamsSet1();

		/* choose exponential map, parallel translation and extrinsic approach and no householder reflections
		Even though the Householder reflections are not used, the locking condition is satisfied.*/
		void ChooseSphereParamsSet2();

		/* choose qf, parallel translation and extrinsic approach and no householder reflections
		The locking conidition is not satisfied*/
		void ChooseSphereParamsSet3();

		/* choose qf, parallel translation and extrinsic approach and no householder reflections
		Beta \neq 1 is used and the locking conidition is satisfied*/
		void ChooseSphereParamsSet4();

		/* choose proximal mapping, parallelization and intrinsic approach and no householder reflections
		the locking conidition is not satisfied. This is used to minimize f(X) + ||X||_1, where f(X) is
		the original cost function defined on a manifold*/
		void ChooseSphereParamsSet5();

		/*Beside the exponential mapping of the sphere, the retractions defined in Stiefel.h also can be used.*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*Beside the exponential mapping of the sphere, the retractions defined in Stiefel.h also can be used.*/
		virtual void Retraction(Variable *x, Vector *etax, Variable *result, double instepsize) const;

		/*Beside the cotangent vector of exponential mapping of the sphere, the retractions defined in Stiefel.h also can be used.*/
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Beside the vector transport by differeitiated the exponential mapping, the differentiated retraction defined in Stiefel.h
		also can be used.*/
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Beside the parallel translation of the unit sphere, the vector transport defined in Stiefel.h also can be used.*/
		virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*Beside the inverse parallel translation of the unit sphere, the inverse vector transport defined in Stiefel.h also can be used.*/
		virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*If parallel translation is used, then call member function "ExpHInvTran",
		otherwise, call function in Stiefel::HInvTran*/
		virtual void HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*If parallel translation is used, then call member function "ExpTranH",
		otherwise, call function in Stiefel::TranH*/
		virtual void TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*If parallel translation is used, then call member function "ExpTranHInvTran",
		otherwise, call function in Stiefel::TranHInvTran*/
		virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);
	private:
		/*Define the proximal update for problems f(x) + ||X||_1.
		The update is result = prox_instepsize(x - etax).*/
		virtual void ProxRetraction(Variable *x, Vector *etax, Variable *result, double instepsize) const;

		/* exponential mapping using extrinsic approach*/
		virtual void ExpRetraction(Variable *x, Vector *etax, Variable *result) const;

		/* the cotangent vector using exponential mapping and extrinsic approach*/
		virtual void ExpcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/* the vector transport by differentiated exponential mapping using extrinsic representation*/
		virtual void ExpDiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/* parallel translation using extrinsic representation*/
		virtual void ExpVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/* inverse parallel translation using extrinsic representation */
		virtual void ExpInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Compute result = H(:, start : end) * \mathcal{T}^{-1}, where H(:, start : end) denotes the matrix formed by columns from "start" to "end".
		\mathcal{T}^{-1} is the inverse parallel translation*/
		virtual void ExpHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H(start : end, :), where H(start : end, :) denotes the matrix formed by cows from "start" to "end".
		\mathcal{T} is the parallel translation*/
		virtual void ExpTranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*Compute result = \mathcal{T} * H * \mathcal{T}^{-1}. \mathcal{T} is the parallel translation*/
		virtual void ExpTranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of SPHERE_H