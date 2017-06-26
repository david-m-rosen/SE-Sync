/*
This file defines the abstract base class for all the solvers
It defines the common properties and features of all the solvers

Solvers

---- WH
*/

#ifndef SOLVERS_H
#define SOLVERS_H

//#include <cmath>
#include "Manifolds/Manifold.h"
#include "Others/def.h"
#include "Problems/Problem.h"
#include <ctime>
#include <iomanip>
#include <iostream>

/*Define the namespace*/
namespace ROPTLIB {

///*Compute min_{y in convex hull of gfs and prefgs are tangent vectors at the
/// tangent space at x} ||y||
// It is defined in SphereConvexHull.h and SphereConvexHull.cpp */
// extern double MinNormConHull(const Manifold *Mani, Variable *x, Vector **Ys,
// integer LYs, Vector *Soln);

/*Compute min_{y in convex hull of gfs and prefgs are tangent vectors at the
tangent space at x} ||y||
It is defined in MinPNormConHull.h and MinPNormConHull.cpp */
extern double MinPNormConHull(const Manifold *Mani, Variable *x, Vector **Ys,
                              integer LYs, Vector *Soln, double *YtY,
                              integer inc);

/*The algorithm is stopped when a value (specified by ther parameter) is less
than the "Tolerance" (a member variable)
The value should be assigned to the member variable: "Stop_Criterion" and the
applicable values are
FUN_REL: fabs((f_k - f_{k+1}) / f_k)
GRAD_F: \|gf_k\|
GRAD_F_0: \|gf_k\| / \|gf_0\|
PSSUBGRAD: min_{Y \in convex hull} conex hull(gf_k, gf_{k-1}, ... )*/
enum StopCrit { FUN_REL, GRAD_F, GRAD_F_0, PSSUBGRAD, STOPCRITLENGTH };

/*Specify what information will be output in the algorithm.
The value should be assigned to the member variable: "Debug",
NOOUTPUT: no output
FINALRESULT: final results are outputted
ITERRESULT: Output information every "OutputGap" iterations, "OutputGap" is a
member variable
DETAILED: Output more than necessary information. Developers can put debug
information for this mode.
The details of output information can be found in Appendix B of the User Manual.
*/
enum DEBUGINFO { NOOUTPUT, FINALRESULT, ITERRESULT, DETAILED, DEBUGLENGTH };

class Solvers {
public:
  /*Run the algorithm. In this class, this function only initialize debug
     information and output the name of algorithm.
          This function has been overloaded for all the algorithms*/
  virtual void Run();

  /*Check whether the general parameters are legal or not.*/
  virtual void CheckParams();

  /*Output all the results after calling the "Run" function.*/
  virtual void OutPutResults(Variable *inx1, double &inf1, double &inngf0,
                             double &inngf, integer &initer, integer &innf,
                             integer &inng, integer &innR, integer &innV,
                             integer &innVp, double &inComTime,
                             double *intimeSeries, double *infunSeries,
                             double *ingradSeries, double *indistSeries,
                             integer &inlengthSeries);

  /*Get the optimizer*/
  inline const Variable *GetXopt(void) const { return x1; };

  /*Get the final cost function value*/
  inline double Getfinalfun(void) const { return f2; };

  /// Added by dmrosen 26-6-2017
  inline double GetPreviousIterateVal() const { return f1; };

  /*Get the norm of the final gradient*/
  inline double Getnormgf(void) const { return ngf; };

  /*Get the norm of the gradient at final iterate over the norm of the gradient
   * at initiate iterate*/
  inline double Getnormgfgf0(void) const { return ngf / ngf0; };

  /*Geth the computational wall time of the algorithm*/
  inline double GetComTime(void) const { return ComTime; };

  /*Get the number of function evaluations*/
  inline integer Getnf(void) const { return nf; };

  /*Get the number of the gradient evaluations*/
  inline integer Getng(void) const { return ng; };

  /*Get the number of retraction evaluations*/
  inline integer GetnR(void) const { return nR; };

  /*The first time to evaluate the action of vector transport
  \mathcal{T}_{\eta_x} is usually
  more expensive than the other times to evaluate the action of vector transport
  \mathcal{T}_{\eta_x}.
  The next two functions are used to
  get the number of action of vector transport (first time)
  get the number of action of vector transport (the other times)
  respectively*/
  inline integer GetnV(void) const { return nV; };
  inline integer GetnVp(void) const { return nVp; };

  /*Get the number of action of Hessian*/
  inline integer GetnH(void) const { return nH; };

  /*Get the number of iterations*/
  inline integer GetIter(void) const { return iter; };

  /*timeSeries, funSeries and gradSeries are three arrays to store the
  computational time, function value and the norm of gradient
  after each iteration
  The next four functions are used to
  get the number of those three series
  get the computation time series
  get the function values series
  get the norm of gradient series
  get the dists from soln to the iterates
  respectively. */
  inline integer GetlengthSeries(void) const { return lengthSeries; };
  inline double *GettimeSeries(void) const { return timeSeries; };
  inline double *GetfunSeries(void) const { return funSeries; };
  inline double *GetgradSeries(void) const { return gradSeries; };
  inline double *GetdistSeries(void) const { return distSeries; };

  /*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e.,
  std::map<std::string, double> .
  This function is used to set the parameters by the mapping*/
  virtual void SetParams(PARAMSMAP params);

  /*Beside the three stopping criterion specified by the member variable
  "Stop_Criterion",
  user also can define a stopping criterion by assigning the following function
  pointer.
  The code always run this function pointer first if it is not a null pointer.
  */
  bool (*StopPtr)(Variable *x, Vector *gf, double f, double ngf, double ngf0,
                  const Problem *prob, const Solvers *solver);

  /*Destructor. It is a pure virtual function*/
  virtual ~Solvers(void) = 0;

  /*member variable: specify the stopping criterion. The applicable values are
  defined in enumerate "StopCrit"
  Default: GRAD_F_0*/
  StopCrit Stop_Criterion;

  /*If a value is less than the Tolerance, then stop the algorithm. The value is
  specified by the member variable "Stop_Criterion".
  Default: 10^{-6}*/
  double Tolerance;

  /*This parameter is used in the stopping criterion for optimizing partly
  smooth functions.
  If ||x_{k+1} - x_k|| / (||x_{k+1}|| + 1) is less than Diffx, then store the
  gradient of latter one.
  Compute the minimum norm vector in the convex hull of the gradients stored. If
  the norm of the minimum norm
  vector is less than Tolerance, then algorithm terminates.
  Default: 1e-6*/
  double Diffx;

  /*The number that is more than the dimension of the domain manifold. This is
  for the solvers of optimizing partly smooth functions.
  In those algorithms, one needs to compute a minimum length in a convex hull of
  a few vectors. The number of those vectors
  are greater than the dimenion of the domain manifold. This parameter specifies
  the number greater than the dimension.
  Default: 3*/
  integer NumExtraGF;

  /*Upper bound for running the algorithm. If excuable time of the algorithm is
  more the the TimeBound, then the algorithm
  is stopped no matther whether a stopping criterion is satisfied or not.
  Default: 60 * 60 * 24 * 365 (one year);*/
  double TimeBound;

  /*Maximum number of iteration. The algorithm will stop if the number of
  iteration reaches Max_Iteration no matter
  whether a stopping criterion is satisfied or not*/
  integer Max_Iteration;

  /*Minimum number of iteration. The algorithm will keep running if the number
  of iteration is less than Min_Iteration
  no matther whether a stopping criterion is satisfied or not*/
  integer Min_Iteration;

  /*If the Debug is larger than FINALRESULT, then the information will be
   * outputted every "OutputGap" iterations*/
  integer OutputGap;

  /*Specified the output information of the algorithm. The applicable values are
   * given in the enumerate "DEBUGINFO".*/
  DEBUGINFO Debug;

protected:
  /*Initialize the solvers by calling the "SetProbX" and "SetDefultParams"
     functions.
          INPUT:	prob is the problem which defines the cost function,
     gradient and possible the action of Hessian
          and specifies the manifold of domain.
          initialx is the initial iterate.
          insoln is the true solution. It is not required and only used for
     research.*/
  virtual void Initialization(const Problem *prob, const Variable *initialx,
                              const Variable *insoln);

  /*Initialize the type of iterates x1, x2 and tangent vectors gf1, gf2 and
     obtian the problem and manifold information
          INPUT:	prob is the problem which defines the cost function,
     gradient and possible the action of Hessian
          and specifies the manifold of domain.
          initialx is the initial iterate.
          insoln is the true solution. It is not required and only used for
     research*/
  virtual void SetProbX(const Problem *prob, const Variable *initialx,
                        const Variable *insoln);

  /*Setting parameters (member variables) to be default values */
  virtual void SetDefaultParams(void);

  /*When one iteration, some algorithms need to update some information. For
     example,
          quasi-Newton methods need to update the Hessian approximation and
     nonlinear conjugate gradient
          needs to update the search direction. They are done in the following
     function*/
  virtual void UpdateData(void) = 0;

  /*Check whether a stopping criterion is satisfied or not*/
  virtual bool IsStopped(void);

  /*Print general information, which is not specific to an algorithm.*/
  virtual void PrintGenInfo(void);

  /*Print information specific to an algorithm*/
  virtual void PrintInfo(void);

  // algorithm-related variables:
  Variable *x1, *x2;   /*x1: current iterate, x2: next iterate*/
  Variable *soln;      /*soln is supposed to be the true solution. If DEBUG >=
                          ITERRESULT and soln != nullptr, then distSeries outputs
                          dist(x_i, soln).*/
  Vector *gf1, *gf2;   /*gf1: gradient at x1, gf2: gradient at x2*/
  double f1, f2;       /*f1: function value at x1, f2: function value at x2*/
  double ngf0, ngf;    /*ngf0: the norm of gradient at the initial iterate, ngf:
                          the norm of the gradient at the final iterate*/
  Vector *eta1, *eta2; /*In Line search-based methods, eta1 is the search
                          direction. eta2 is stepsize * eta1.*/
  /*In trust region-based methods, eta1 is the initial guess for solving the
  local model.
  eta2 is the result produced by solving the local model.*/
  Vector *zeta; // temp storage

  // the next six parameters are for stopping criterion of optimizing partly
  // smooth cost functions
  Vector **gfs;             /*The gradients of previous few steps*/
  integer Lengthgfs;        /*The maximum length of gfs*/
  integer Currentlengthgfs; /*the current of gfs*/
  integer idxgfs;           /*current idx of gfs*/
  double nsubgf; /*the norm of the minimum gradient in the convex hall*/
  integer
      subprobtimes; /*the number of runs for solving the sub convex problem*/

  /*points in a neighborhood of current iterate. It is used in RGS.*/
  Variable **Xs;

  // Input parameters and functions
  const Manifold *Mani; /*The manifold on which the cost function is*/
  const Problem *Prob;  /*The problem which defines the cost function, gradient
                           and probably action of Hessian*/

  // For debug information
  integer iter;            /*number of iterations*/
  unsigned long starttime; /*the start time of running the algorithm*/
  double ComTime;          /*the computational time*/
  integer nf, ng, nR, nV, nVp,
      nH; /*number of function evaluations
                                          number of gradient evaluations
                                          number of retraction evaluations
                                          number of vector transport (See
             GetnV(void) for details)
                                          number of vector transport (See
             GetnVp(void) for details)
                                          number of action of Hessian*/
  double *timeSeries, *funSeries, *gradSeries,
      *distSeries;      /*four arrays to store the computational time, function
                           values, norm of gradient
                                                       after each iteration, and
                           dist(soln, x_i)*/
  integer lengthSeries; /*the length of above four arrays, i.e., the length of
                           timeSeries, funSeries, gradSeries, distSeries.*/
  std::string SolverName; /*The name of the solver. This is assigned in the
                             constructor function of each derived class*/

  /*new memory for the double array Vs, type Vector, with length l*/
  void NewVectors(Vector **&Vs, integer l);

  /*delete memory for the double array Vs, type Vector, with length l*/
  void DeleteVectors(Vector **&Vs, integer l);

  /*new memory for the double array Xs, type Variable, with length l*/
  void NewVariables(Vector **&Xs, integer l);

  /*delete memory for the double array Xs, type Variable, with length l*/
  void DeleteVariables(Vector **&Xs, integer l);
};

}; /*end of ROPTLIB namespace*/

#endif // end of SOLVER_H
