/*
Driver for the problem defines in [HGSA2015]
	[HGSA2015]: Wen Huang and K. A. Gallivan and Anuj Srivastava and P.-A. Absil, Riemannian Optimization for Registration of Curves in Elastic Shape Analysis
		Journal of Mathematical Imaging and Vision, doi = "10.1007/s10851-015-0606-8", 2015
Find best reparameterization (rotation) for functions, open and closed curves in R^d
\int_D \|Oq1(t) - q2(\int_0^t l^4(s) ds + m mod 1) l^2(t)\|_2^2 dt + w B(l)
the domain D = [0, 1] for open curves, D is a circle for closed curves.

This driver does not only find a local minimizer. It uses multiple initial points to obtain
acceptable global optimum. The idea of computing initial points is given in [HGSA2015].

---- WH
*/

#ifndef DRIVERELASTICCURVESRO_H
#define DRIVERELASTICCURVESRO_H


#include <iostream>
#include "Others/randgen.h"
#include "Manifolds/Manifold.h"
#include "Problems/Problem.h"
#include "Solvers/SolversLS.h"
#include <ctime>

#include "Manifolds/Euclidean/EucVariable.h"
#include "Manifolds/Euclidean/EucVector.h"

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

#include "Problems/ElasticCurvesRO/ElasticCurvesRO.h"
#include "Others/MyMatrix.h"

#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*The used slopes in Dynamic Programming*/
#define NNBRS 23
	const integer Nbrs[NNBRS][2] = {
		{ 1, 1 },
		{ 1, 2 },
		{ 2, 1 },
		{ 2, 3 },
		{ 3, 2 },
		{ 1, 3 },
		{ 3, 1 },
		{ 1, 4 },
		{ 3, 4 },
		{ 4, 3 },
		{ 4, 1 },
		{ 1, 5 },
		{ 2, 5 },
		{ 3, 5 },
		{ 4, 5 },
		{ 5, 4 },
		{ 5, 3 },
		{ 5, 2 },
		{ 5, 1 },
		{ 1, 6 },
		{ 5, 6 },
		{ 6, 5 },
		{ 6, 1 }
	};

	/*Indicates how many slopes used in Dynamic Programming
	NUMBIG: use all 23 choices.
	NUMSMALL: use 11 choices
	*/
	enum SLOPESTYPE{ NUMBIG, NUMSMALL };

	/*The entrance of the driver
	C1, C2 are two curves in R^d. The curves are represented by n points.
	w is the coefficient used in the barrier function.
	rotated: whether rotation is considered or not.
	isclosed: whether the curve is closed or not.
	onlyDP: Dynamic Programming is used without using Riemannian optimization to improve the solution.
	skipm: the interval between two break points is at least skipm.
	*/
	void DriverElasticCurvesRO(double *C1, double *C2, integer d, integer n, double w, bool rotated, bool isclosed,
		bool onlyDP, integer skipm, std::string solverstr, integer autoselectC, ProductElement *Xopt, bool &swap, double *fopts, 
		double *comtime, integer &Nsout, integer &numinitialx, double *optQ1 = nullptr, double *optQ2 = nullptr);

	/*The dynamic programming*/
	double DynamicProgramming(const double *p_q1, const double *p_q2, integer d, integer N, double *gamma, bool isclosed, SLOPESTYPE Nbrstype);

	/*Centerize the curve C, i.e., make the mean of the curve in origin*/
	void CenterC(double *C, integer d, integer n);

	/*Make the two norm of C equal to 1*/
	void NormalizedC(double *C, integer d, integer n);

	/*Compute the total change of the angle along the curve C*/
	double ComputeTotalAngle(const double *C, integer d, integer n);

	/*Find the initial break points and the Ns for the coarse dynamic programming*/
	void FindInitialBreaksAndNs(const double *C, integer d, integer n, integer minSkip, double thresholdsmall,
		integer rand_shift, integer *p_ms, integer &Lms, integer &Ns);

	/*Compute the q (SRV) function from the curve*/
	void CurveToQ(const double *C, integer d, integer n, double *q, bool isclosed);

	/*Compute the curve from its q (SRV) function*/
	void QToCurve(const double *Q, integer d, integer n, double *C, bool isclosed);

	/*Shift the break point by m points.*/
	void ShiftC(const double *C, integer d, integer n, double *Cshift, integer m);

	/*Find the best rotation between q1 and q2, i.e., \min_O \|Oq1 - q2\|_{L^2}*/
	void FindBestRotation(const double *q1, const double *q2, integer d, integer n, double *O);

	/*Resample the curve using cubic spline interpolation and use ns points to represent the curve*/
	void GetCurveSmall(const double *C, double *Cs, integer d, integer n, integer ns, bool isclosed);

	/*Compute the inverse of the diffeomorphism gamma: [0, 1] --> [0, 1]*/
	void GammaInverse(const double *DPgam, integer n, double *DPgamI);

	/*Resample the gamma function using linear interploation, and using n points to represent the gamma function*/
	void ReSampleGamma(const double *DPgams, integer ns, double *DPgam, integer n);

	/*Compute the gradient using finite difference (Center)*/
	void Gradient(const double *DPgam, integer n, double h, double *grad);

	/*Compute the gradient using finite difference (Center)*/
	void GradientPeriod(const double *DPgam, integer n, double h, double *grad);

}; /*end of ROPTLIB namespace*/

#endif // end of DRIVERELASTICCURVESRO_H
