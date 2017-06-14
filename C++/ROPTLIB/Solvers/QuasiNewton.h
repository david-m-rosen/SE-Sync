/*
This file defines the abstract base class for the Hessian approximation
operations in quasi-Newton methods.

Solvers --> QuasiNewton

---- WH
*/

#ifndef QUASINEWTON_H
#define QUASINEWTON_H

//#include <cmath>
#include <iostream>
#include <list>
#include <iomanip>
#include <ctime>
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/Solvers.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	class QuasiNewton:public Solvers{
	public:
		/*specify whether the cost function is convex or not.
		If yes, then the initial Hessian approximation is a scalar times identity, where the scalar is to
		measure the magnitude of eigenvalues, otherwise, it is identity.
		Default: false*/
		bool isconvex;

		/*The same as \epsilon in [LF01, (3.2)]
		[LF01]: D.-H. Li and M. Fukushima. On the global convergence of the BFGS method for nonconvex unconstrained optimization problems.
		SIAM Journal on Optimization, 11(4):1054?064, 2001
		Default: 10^{-4}*/
		double nu;

		/*The same as \alpha in [LF01, (3.2)]
		[LF01]: D.-H. Li and M. Fukushima. On the global convergence of the BFGS method for nonconvex unconstrained optimization problems.
		SIAM Journal on Optimization, 11(4):1054?064, 2001
		Default: 1*/
		double mu;

		/*the number of pairs of s and y in Limited-memory version of quasi-Newton methods;  The same as \ell in [HGA2015, Algorithm 2]
		[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685, 2015.
		Default: 4*/
		integer LengthSY;

		/*the number of previous bb1 stepsize. Used in ABB_min stepsize. See details in [SRTZ2017]
		[SRTZ2017]: On the steplength selection in gradient methods for unconstrained optimization.
		stepsize * id can be used as the initial Hessian approximation in limite-memory quasi-Newton methods
		Default: 0*/
		integer Num_pre_BB;

		/*ratio for step size selection. It is the same as \tau in [Algorithm2, SRTZ2017]
		stepsize * id can be used as the initial Hessian approximation in limite-memory quasi-Newton methods
		[SRTZ2017]: On the steplength selection in gradient methods for unconstrained optimization.
		Default: 0, which defines ss/sy stepsize, if it is 1 and Num_pre_BB1 is 0, then it defines sy/yy stepsize.*/
		double BBratio;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Compute result = H v in RBroydenFamily*/
		virtual void HvRBroydenFamily(Vector *v, Vector *result);

		/*Update the Hessian approximation for RBroydenFamily if necessary*/
		virtual void UpdateDataRBroydenFamily(void);

		/*The phi coefficient in [HGA2015, (2.3)]
		[HGA2015]: Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685, 2015
		*/
		double Phi(Variable *x2, Vector *y, Vector *s, LinearOPE *tildeH, double inpsy, double yHy, Vector *u);

		/*Compute result = H v in RBFGS*/
		virtual void HvRBFGS(Vector *v, Vector *result);

		/*Update the Hessian approximation for RBFGS if necessary*/
		virtual void UpdateDataRBFGS(void);

		/*Compute result = H v in LRBFGS*/
		virtual void HvLRBFGS(Vector *v, Vector *result);

		/*Update the Hessian approximation for LRBFGS if necessary*/
		virtual void UpdateDataLRBFGS(void);

		/*Compute result = H v in RWRBFGS*/
		virtual void HvRWRBFGS(Vector *v, Vector *result);

		/*Update the Hessian approximation for RWRBFGS if necessary*/
		virtual void UpdateDataRWRBFGS(void);

		/*Compute result = H v in RTRSR1*/
		virtual void HvRTRSR1(Vector *v, Vector *result);

		/*Update the Hessian approximation for RTRSR1 if necessary*/
		virtual void UpdateDataRTRSR1(void);

		/*Compute result = H v in LRTRSR1*/
		virtual void HvLRTRSR1(Vector *v, Vector *result);

		/*Update the Hessian approximation for LRTRSR1 if necessary*/
		virtual void UpdateDataLRTRSR1(void);

		/*Compute result = H v in RBFGS for L-continuous functions*/
		virtual void HvRBFGSSub(Vector *v, Vector *result);

		/*Update the Hessian approximation for RBFGS for L-continuous functions if necessary*/
		virtual void UpdateDataRBFGSSub(void);

		/*Compute result = H v in LRBFGS for L-continuous functions*/
		virtual void HvLRBFGSSub(Vector *v, Vector *result);

		/*Update the Hessian approximation for LRBFGS for L-continuous functions if necessary*/
		virtual void UpdateDataLRBFGSSub(void);

		/*===================RBFGS subgradient Lipschitz continous====================*/
		/*lower bound of the smallest eigenvalue of the hessian approximation
		Default: 1e-7*/
		double lambdaLower;
		/*upper bound of the largest eigenvalue of the hessian approximation
		Default: 1e7*/
		double lambdaUpper;
	protected:

		/*===================LRBFGS====================*/
		/*initial Hessian approximation in limited-memory BFGS method. It is a scalar times identity.*/
		virtual double InitialHessian(double inpss, double inpsy, double inpyy);
		std::list<double> pre_BBs; /* Store a few computed BB stepsize ss/sy for initial Hessian approximation using adaptive BB min (ABB_min) idea*/

		/*===================RBFGS, RBroydenfamily, RTRSR1====================*/
		bool isupdated; /*Mark whether the (inverse) Hessian approximation is updated*/
		double betay, phic, inpsy, inpss, inpyy;  /*betay: \|\xi\| / \|\mathcal{T}_{R_\xi} \xi\| in the locking condition;
												phic: the coefficient (1-phic) BFGS + phi DFP in Broyden family method
												inpsy: g(s, y); inpss: g(s, s); inpyy: g(y, y); */

		Vector *s, *y, *u;/*the s, y, and u of current step*/
		LinearOPE *H, *tildeH; /*The inverse Hessian approximations for current and next iterations respectively*/
		LinearOPE *B, *tildeB; /*The Hessian approximations for current and next iterations respectively*/

		/*===================LRBFGS LRTRSR1====================*/
		Vector **S, **Y; /*The stored pairs of s and y*/
		double *RHO; /*the sequence of 1 / g(s_k, y_k), where Hessian approximation k-th iteration is updated*/
		double rho, gamma; /*rho: 1/g(s, y) at current iteration, gamma: g(s, y) / g(y, y) for LRBFGS and gamma: g(y, y) / g(s, y) for RTRSR1*/
		integer Currentlength; /*The current length of array S, Y and RHO*/
		integer beginidx; /*The starting index of S, Y and RHO at current iteration*/

		/*===================LRTRSR1====================*/
		bool ischangedSandY; /*Mark whether S and Y is updated.*/
		Vector **YMGS; /*The stored pairs of s and y, and also Y - gamma S*/
		double *SS, *SY, *PMGQ; /*SS is S^\flat S which is the matrix Q defined in [HAG2014, (46)],
								SY is the matrix P defined in [HAG2014, (46)]
								PMGQ is the P - gamma Q matrix defined in [HAG2014, (46)] */
		integer *P;	/*the pemuation matrix when computing the LU decomposition for the matrix PMGQ*/

	};
}; /*end of ROPTLIB namespace*/

#endif // end of QUASINEWTON_H
