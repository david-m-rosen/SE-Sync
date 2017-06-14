/*
This file defines the class of the line search algorithm for locally lipschitz functions on Riemannian manifolds

Solvers --> QuasiNewton --> SolversLS --> SolversLSLPSub

---- WH
*/

#ifndef SOLVERSLPSUB_H
#define SOLVERSLPSUB_H

#include "Solvers/SolversLS.h"
#include "Others/MinPNormConHull.h"
#include "Manifolds/Sphere/Sphere.h"
#include "Solvers/LRBFGS.h"
#include "Solvers/RTRNewton.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class SolversLSLPSub : public SolversLS{
	public:
		/*Run the algorithm. This function gives the framework for the linesearch method*/
		virtual void Run(void);

		/*Destructor. Delete the vectors and Hessian approximation used in RBFGSLPSub, i.e., s and y, H and tildeH*/
		virtual ~SolversLSLPSub(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Check whether the parameters about RBFGSLPSub are legal or not.*/
		virtual void CheckParams(void);

		/*Initialize the type of iterates x1, x2 and tangent vectors gf and obtian the problem and manifold information
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*Check whether a stopping criterion is satisfied or not*/
		virtual bool IsStopped(void);

		/*Epsilon-subdifferential
		Default: 1*/
		double Eps;
		/*Reduce the Epsilon by Theta_eps each time, i.e., Eps *= Theta_eps. 
		Default: 0.01*/
		double Theta_eps;
		/*The minimum value of Eps. Too small epsilon have non-negligible numerical error.
		Default: 1e-6*/
		double Min_Eps;

		/*The bound on the square of P-norm of gradient 
		Default: 1*/
		double Del;
		/*Reduce Del by Theta_del each time.
		Default: 0.01*/
		double Theta_del;

	protected:
		double MinPnormWv(void);
		void Increasing(double neta1, double Pngfsq, double hb);

		/*Compute the search direction */
		void GetSearchDir(void);

		/*A function pointer to an action of the Hessian*/
		void (QuasiNewton::*Hv)(Vector *v, Vector *result);

		/*Algorithm-related variables*/
		Vector *gf; /*the subgradient with minimum P norm in eps-subgradient of f*/

		double *Sub_soln; /*The solution of the convex hull optimization problem.
						  It can use as the initial iterate for next run.*/

		//double Pngfsq;
		//double neta1;
		//double ht, hb;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of SOLVERSLPSUB_H
