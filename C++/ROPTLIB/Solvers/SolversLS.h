/*
This file defines the abstract base class for all the linesearch-based solvers
It defines the common properties and features of all the linesearch-based solvers

Solvers --> QuasiNewton --> SolversLS
							
---- WH
*/

#ifndef SOLVERLS_H
#define SOLVERLS_H

//#include <cmath>
#include <iostream>
#include <list>
#include <ctime>
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/Solvers.h"
#include "Solvers/QuasiNewton.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/* Linesearch algorithms. It should be assigned to the member variable "LineSearch_LS".
	ARMIJO: The Armijo-Goldstein condition. [DS83 Algorithm A6.3.1]
	WOLFE:  The weak Wolfe condition. [DS83 Algorithm A6.3.1mod]
	STRONGWOLFE: The strong Wolfe condition. [NW06 Algorithm 3.5]
	EXACT: The exact line search based on scalar quasi-Newton method
	WOLFELP: The weak Wolfe condition for lipschitz continuous functions
	INPUTFUN: For this option, users can specify their own line search algorithm by assigning the function pointer "LinesearchInput".
	[DS83]: J. E. Dennis and R. B. Schnabel. Numerical methods for unconstrained optimization and nonlinear equations. Springer, New Jersey, 1983
	[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006
	*/
	enum LSAlgo{ ARMIJO, WOLFE, STRONGWOLFE, EXACT, WOLFELP, INPUTFUN, LSALGOLENGTH };

	/* Linesearch status. It is an output argument and users don't need to assign this enumerate to any member variable.
	NOCURVATURE: the second Wolfe condition is not satisfied
	MINSTEPSIZE: line search algorithm reaches the minimum stepsize
	MAXSTEPSIZE: line search algorithm reaches the maximum stepsize
	NONEXACT: exact line search algorithm does not find a point satisfying the inner stopping criterion
	SUCCESS: line search algorithm succeeds in finding a point satisfying the line search condition
	*/
	enum LSstatusSet{ NOCURVATURE, MINSTEPSIZE, MAXSTEPSIZE, NONEXACT, LSERROR, SUCCESS, LSSTATUSSETLENGTH };

	/*Initial step size in line search algorithm.
	ONESTEP: t0 = one 
	BBSTEP: t0 = g(s, s) / g(s, y), s is the difference of consecutive iterates and y is the difference of the
			gradients at consecutie iterates.
	QUADINT: t0 = [(3.60), NW06]
	QUADINTMOD: t0 = [page 60, NW06]
	[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006
	*/
	enum InitStepsizeSet{ ONESTEP, BBSTEP, QUADINT, QUADINTMOD, EXTRBBSTEP, INITSTEPSIZESETLENGTH };

	class SolversLS : public QuasiNewton{
	public:
		/*Run the algorithm. This function gives the framework for all the linesearch based methods*/
		virtual void Run(void);

		/*Check whether the parameters about linesearch algorithms are legal or not.*/
		virtual void CheckParams(void);

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Beside the four line search algorithms provided in this library and specified by the member variable "LineSearch_LS",
		user also can define a line search algorithm by assigning the following function pointer.
		User needs to assign LineSearch_LS to be INPUTFUN to call this function. */
		double(*LinesearchInput)(integer iter, Variable *x1, Vector *eta1, double initialstepsize, double initialslope, const Problem *prob, const Solvers *solver);

		/* ===============public parameters below================= */

		/*Line search algorithm. The applicable values are in the enumerate LSAlgo
		Default: ARMIJO */
		LSAlgo LineSearch_LS;

		/*If IsPureLSInput is false, then the step size from the user-written line search method, will
		be used as the initial stepsize for the Armijo linesearch. Otherwise, the step size from the user-written line search
		method is used as the final step size.
		Default: false*/
		bool IsPureLSInput;

		/*the coefficient of the Wolfe first condition
		Default: 0.0001 */
		double LS_alpha;

		/*the coefficient of the Wolfe second condition
		Default: 0.999*/
		double LS_beta;

		/*the minimum stepsize allowed in the linesearch algorithm
		Default: machine eps*/
		double Minstepsize;

		/*the maximum stepsize allowed in the linesearch algorithm
		Default: 1000 */
		double Maxstepsize;

		/*The coefficient in the Armijo-Goldstein condition
		Default: 0.1 */
		double LS_ratio1;

		/*The coefficient in the Armijo-Goldstein condition
		Default: 0.9 */
		double LS_ratio2;

		/*Initial stepsize at the first iteration
		Default: 1*/
		double Initstepsize;

		/*If the norm of current gradient over the norm of the initial gradient is less than Accuracy,
		then the step size is fixed to be the member variable "Finalstepsize".
		Defaut: 0*/
		double Accuracy;

		/*When the iterate is close to the minimizer (see the annotation of Accuracy), then fixed
		the stepsize to the Finalstepsize if Finalstepsize > 0. If Finalstepsize <= 0, then the proposed
		initial stepsize is used as the accepted stepsize.
		Default: 1*/
		double Finalstepsize;

		/* the number of computed functions values. This is used in the nonmonotonic linesearch which
		uses max_{1 \leq i \leq num_pre_funs} (f_{k + 1 - i})
		Default: 0, which defines the standard Armijo line search.*/
		integer Num_pre_funs;

		/*Line search algorithm. The applicable values are in the enumerate InitStepsizeSet
		Default: QUADINTMOD (It may alter based on the derived algorithm class) */
		InitStepsizeSet InitSteptype;
	protected:

		/*initial Hessian approximation in limited-memory BFGS method. It is a scalar times identity.*/
		virtual double InitialHessian(double inpss, double inpsy, double inpyy);

		/*Choose what line search algorithm is used*/
		virtual void ChooseLinesearch(void);

		/*Print general information of linesearch based algorithms, which is not specific to an algorithm.*/
		virtual void PrintGenInfo(void);

		/*Delete objects that are used in this class*/
		virtual ~SolversLS(void);

		/*Compute the search direction. It is a pure virtual function.*/
		virtual void GetSearchDir(void) = 0; // required to be overload in derived class if the derived class is not abstract

		/*Compute the initial stepsize using [NW06, Page 60]
			[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006	*/
		virtual void InitialStepSize(void);

		/*Call Solvers::SetProbX function and set up the temporary objects for linesearch algorithm.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams(void);

		/*The Armijo-Goldstein condition.[DS83 Algorithm A6.3.1] combined with nonmontone line search
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void LinesearchArmijo(void);

		/*The weak Wolfe condition [DS83 Algorithm A6.3.1mod]
		[DS83] : J.E.Dennis and R.B.Schnabel.Numerical methods for unconstrained optimization and nonlinear equations.Springer, New Jersey, 1983*/
		virtual void LinesearchWolfe(void);

		/*The strong Wolfe condition [NW06 Algorithm 3.5]
		[NW06] : J.Nocedal and S.J.Wright.Numerical optimization.Springer, second edition, 2006*/
		virtual void LinesearchStrongWolfe(void);

		/*The exact line search based on scalar quasi-Newton method*/
		virtual void LinesearchExact(void);

		/*The weak Wolfe condition for Lipschitz functions*/
		virtual void LinesearchWolfeLipschitz(void);

		/*Evaluate the cost function h(stepsize) = f(R_{x_1}(stepsize * eta1))*/
		virtual double h(void);

		/*Evaluate the derivative of cost function h, i.e., h'(stepsize) = \frac{d}{d stepsize} f(R_{x_1}(stepsize * eta1))*/
		virtual double dh(void);

		/*When one iteration, some algorithms need to update some information. For example,
		quasi-Newton methods need to update the Hessian approximation and nonlinear conjugate gradient
		needs to update the search direction. They are done in the following function*/
		virtual void UpdateData(void);

		/*The "Run" function call this function pointer and the function pointer points to one of the member linesearch functions:
			void LinesearchArmijo();
			void LinesearchWolfe();
			void LinesearchStrongWolfe();
			void LinesearchExact(); */
		void (SolversLS::*Linesearch)(void);

		// parameters
		double initiallength;	/*The initial stepsize at an iteration*/
		double stepsize;		/*The step size*/
		double initialslope, newslope;	/*The slopes for the scalar function h(t) = f(R_{x_1}(t * eta1)) at 0 and at the accepted stepsize*/
		std::list<double> pre_funs; /* Store a few computed function values for nonmonotonic line search*/
		LSstatusSet LSstatus;	/*The line search status produced by the linesearch algorithm*/
		std::string *LSstatusSetnames;	/*This string array is to store the line search status names*/
		Vector *exeta1, *exeta2; /*If the intrinsic representation is used, then exeta1 and exeta2 are used to store the extrinsic representations of eta1 and eta2 respectively*/
	private:
		/*The function used in the strong Wolfe condition. See Algorithm 3.6 in [NW06]
		[NW06] : J.Nocedal and S.J.Wright.Numerical optimization.Springer, second edition, 2006 */
		void Zoom(double x1, double fx1, double slopex1, double x2, double fx2);
		void ZoomLP(double a, double b);
	};
}; /*end of ROPTLIB namespace*/
#endif // end of SOLVERLS_H
