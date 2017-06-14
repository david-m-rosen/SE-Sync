/*
This file defines the abstract base class for all the trust region-based solvers
It defines the common properties and features of all the trust region-based solvers

Solvers --> QuasiNewton --> SolversTR

---- WH
*/

#ifndef SOLVERSTR_H
#define SOLVERSTR_H

#include "Solvers/QuasiNewton.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/* output status of the truncated conjugate gradient. It is an output argument and users don't need to assign this enumerate to any member variable.
	TR_NEGCURVTURE: Find negative curvature
	TR_EXCREGION: The resulting vector is out of the trust region
	TR_LCON: Teminate when the kappa variable takes effect, which indicate the convergence rate is linear.
	TR_SCON: Teminate when the theta variable takes effect, which indicate the convergence rate is superlinear.
	TR_MAXITER: Teminate when the inner iterations reach the maximum inner iterations specified by the member variable "Max_Inner_Iter"
	*/
	enum tCGstatusSet{ TR_NEGCURVTURE, TR_EXCREGION, TR_LCON, TR_SCON, TR_MAXITER, TCGSTATUSSETLENGTH };

	class SolversTR : public QuasiNewton{
	public:
		/*Run the algorithm. This function gives the framework for all the trust region based methods*/
		virtual void Run(void);

		/*Check whether the parameters about trust region algorithms are legal or not.*/
		virtual void CheckParams(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/* ===============public parameters below================= */

		/*if the difference between the local model and the true function is greater than "Acceptence_Rho",
		then accept the candadite iterate. See c in [Step 5 in Algorithm 1, HAG2014].
		[HAG2014]:  W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method. Mathematical Programming.
		Default: 0.1 */
		double Acceptence_Rho;

		/*if the local model and the true function does not match well, then the radius is shrinked by "Shrinked_tau",
		See [Step 17 in Algorithm 1, HAG2014]
		Default: 0.25*/
		double Shrinked_tau;

		/*if the local model and the true function match good enough, then the radius is magnified by "Magnified_tau",
		See [Step 12 in Algorithm 1, HAG2014]
		Default: 2*/
		double Magnified_tau;

		/*Allowed minimum radius of the trust region.
		Default: machine eps*/
		double minimum_Delta;

		/*Allowed maximum radius of the trust region.
		Default: 1000*/
		double maximum_Delta;

		/*if useRand is true, then set r to be Hess f(x1)[eta1] + grad f(x1), otherwise, set r to be grad f(x1).
		Default: false*/
		bool useRand;

		/*the maximum iterations allowed for solving the local trust region
		Default: 1000*/
		integer Max_Inner_Iter;

		/*the minimum iterations allowed for solving the local trust region
		Default: 0*/
		integer Min_Inner_Iter;

		/*the theta and kappa are used to check whether stopping criterion of solving local model is satisfied or not
		Default values depends on methods. See the User Manual*/
		double theta;
		double kappa;

		/*The initial radius of the trust region
		Default: 1*/
		double initial_Delta;
	protected:
		/*Print general information of trust region based algorithms, which is not specific to an algorithm.*/
		virtual void PrintGenInfo(void);

		/*Delete objects that are used in this class*/
		virtual ~SolversTR(void);

		/*Compute result = H[Eta], where H is the Hessian or the Hessian approximation*/
		virtual void HessianEta(Vector *Eta, Vector *result) = 0; // required to be overloaded in derived class

		/*Compute the initial iterate. The default value is a zero vector.*/
		virtual void InitialVector(void);

		/*Run the truncated conjugate gradient method for the local model*/
		virtual void tCG_TR(void);

		/*Preconditioner for the solving the local model*/
		virtual void PreConditioner(Variable *x, Vector *eta, Vector *result);

		/*Call Solvers::SetProbX function and set up the temporary objects for trust region algorithm.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*When one iteration, some algorithms need to update some information. For example,
		quasi-Newton methods need to update the Hessian approximation. They are done in the following function*/
		virtual void UpdateData(void);

		/*This function is called when the candidate is accepted. It can be used to transport tangent vectors and Hessian approximation
		from the tangent space at x1 to the tangent space at x2.*/
		virtual void Acceptence(void);

		// algorithm-related variables:
		Vector *r, *z, *delta, *Hd; /*Used for solving the local model*/
		double rho;	/*the difference between the local model and the true function*/
		double Delta;	/*the radius of the trust region*/
		integer innerIter;	/*The number of inner iterations for solving the local model.*/
		tCGstatusSet tCGstatus; /*The status of solving the local model*/
		std::string *tCGstatusSetnames;	/*This string array is to store the trust region status names*/
	};
}; /*end of ROPTLIB namespace*/
#endif
