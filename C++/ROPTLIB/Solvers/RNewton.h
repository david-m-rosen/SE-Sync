/*
This file defines the class of the Riemannian Newton method. This code does not follow a particular paper.
The search direction is by applying truncated conjugate gradient to approximately solve the linear
system Hessian[direction] = - gradient.

Solvers --> SolversLS --> RNewton

---- WH
*/

#ifndef RNEWTON_H
#define RNEWTON_H

#include "Solvers/SolversLS.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/* output status of the truncated conjugate gradient. It is an output argument and users don't need to assign this enumerate to any member variable.
	LS_NEGCURVTURE: Find negative curvature
	LS_LCON: Teminate when the kappa variable takes effect, which indicate the convergence rate is linear.
	LS_SCON: Teminate when the theta variable takes effect, which indicate the convergence rate is superlinear.
	LS_MAXITER: Teminate when the inner iterations reach the maximum inner iterations specified by the member variable "Max_Inner_Iter"
	*/
	enum tCGLSstatusSet{ LS_NEGCURVTURE, LS_LCON, LS_SCON, LS_MAXITER, TCGLSSTATUSSETLENGTH };

	class RNewton : public SolversLS{
	public:
		/*The contructor of RNewton method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		RNewton(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Check whether the parameters about RNewton are legal or not.*/
		virtual void CheckParams();

		/*Destructor. Delete temporary Vectors r, z, delta and Hd; Delete the strings of RCGmethods' names*/
		virtual ~RNewton();

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Call Solvers::SetProbX function; initialize temporary vectors; and indicate RNewton need action of Hessian.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

		/* ===============public parameters below================= */

		/*if useRand is true, then set r to be Hess f(x1)[eta1] + grad f(x1), otherwise, set r to be grad f(x1).
		Default: false*/
		bool useRand;

		/*the maximum iterations allowed for solving the local linear system
		Default: 1000*/
		integer Max_Inner_Iter;

		/*the minimum iterations allowed for solving the local linear system
		Default: 0*/
		integer Min_Inner_Iter;

		/*the theta and kappa are used to check whether stopping criterion of solving local linear system is satisfied or not
		Default: theta 0.1, kappa 0.9*/
		double theta;
		double kappa;
	protected:
		/*Compute the search direction by using truncated conjugate gradient to approximately solve Hessian[direction] = - gradient*/
		virtual void GetSearchDir();

		/*Preconditioner for the solving the linear system Hessian[direction] = - gradient*/
		void PreConditioner(Variable *x, Vector *eta, Vector *result);

		/*Compute result = H[Eta], where H is the Hessian*/
		void HessianEta(Vector *Eta, Vector *result) const;

		/*Print information specific to RNewton*/
		virtual void PrintInfo();

		/*Run the truncated conjugate gradient method for the local linear system*/
		void tCG_LS();

		Vector *r, *z, *delta, *Hd; /*temporary vectors used in the trancated conjugate gradient method*/
		integer innerIter; /*the number of iterations of the truncated conjugate gradient method*/
		tCGLSstatusSet tCGLSstatus; /*the output status of the truncated conjugate gradient method*/
		std::string *tCGLSstatusSetnames; /*the output names of the truncated conjugate gradient method*/
	};
}; /*end of ROPTLIB namespace*/
#endif // end of RNEWTON_H
