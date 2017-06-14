/*
This file defines the class of the Riemannian nonlinear conjugate gradient method. This code does not follow a particular paper.
It simply generates all the existing Euclidean nonlinear conjugate gradient methods to the Riemannian setting.

Solvers --> SolversLS --> RCG

---- WH
*/

#ifndef RCG_H
#define RCG_H

#include "Solvers/SolversLS.h"
#include "Others/def.h"
#undef max

/*Define the namespace*/
namespace ROPTLIB{

	/* Riemannian nonlinear conjugate gradient formulas. It should be assigned to the member variable "RCGmethod".
	The Euclidean formulas can be found in e.g., [NW06, Section 5.2].
	[NW06]: J. Nocedal and S. J. Wright. Numerical optimization. Springer, second edition, 2006
	*/
	enum RCGmethods{ FLETCHER_REEVES, POLAK_RIBIERE_MOD, HESTENES_STIEFEL, FR_PR, DAI_YUAN, HAGER_ZHANG, RCGMETHODSLENGTH };

	class RCG : public SolversLS{
	public:
		/*The contructor of RCG method. It calls the function Solvers::Initialization.
		INPUT : prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		RCG(const Problem *prob, const Variable *initialx, const Variable *insoln = nullptr);

		/*Destructor. Delete the strings of RCGmethods' names*/
		virtual ~RCG();

		/*Check whether the parameters about RCG are legal or not.*/
		virtual void CheckParams();

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);

		/*Call Solvers::SetProbX function and indicate RCG does not need action of Hessian.
		INPUT:	prob is the problem which defines the cost function, gradient and possible the action of Hessian
		and specifies the manifold of domain.
		initialx is the initial iterate.
		insoln is the true solution. It is not required and only used for research.*/
		virtual void SetProbX(const Problem *prob, const Variable *initialx, const Variable *insoln);

		/*Setting parameters (member variables) to be default values */
		virtual void SetDefaultParams();

		/* ===============public parameters below================= */

		/*Reset the search direction to be the negative gradient every "ManDim" iterations. Ideally, "ManDim" shoud be the dimension of the domain manifold.*/
		integer ManDim;

		/*Indicate what formula is used in RCG method*/
		RCGmethods RCGmethod;

	protected:
		/*Compute the search direction based on the RCG forumla.
		Reset the search direction to be negative gradient if number of iterations mod "ManDim" is zero and the search direction is not sufficiently
		descent.*/
		virtual void GetSearchDir();

		/*Compute a candadite of the search direction*/
		virtual void UpdateData();

		/*Print information specific to RCG*/
		virtual void PrintInfo();

		/*Strings to store the names of formula used in RCG methods*/
		std::string *RCGmethodSetnames;

		/*sigma is the coefficient in - \grad f(x_{k+1}) + sigma eta_k*/
		double sigma;
	};
}; /*end of ROPTLIB namespace*/
#endif // end of RCG_H
