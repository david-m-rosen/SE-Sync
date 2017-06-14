/*
This file defines the class for the Stiefel manifold \St(p, n) = \{X \in R^{n \times p} | X^T X = I_p\}
It defines the common properties and features of the manifold.

Manifold --> Stiefel

---- WH
*/

#ifndef STIEFEL_H
#define STIEFEL_H

//#include <cmath>

#include "Manifolds/Stiefel/StieVariable.h"
#include "Manifolds/Stiefel/StieVector.h"
#include "Manifolds/Manifold.h"
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	/*Note that not all metrics, retractions and vector transports have been done.*/

	/* Riemannian Metric for the Stiefel manifold:
	Eucldean: g_x(etax, xix) = \trace(etax^T xix);
	Canonical: g_x(etax, xix) = \trace(etax^T (I_n - x x^T / 2) xix); */
	enum StieMetric{ EUCLIDEAN, CANONICAL, STIEMETRICLENGTH };

	/*Retraction for the Stiefel manifold
	QF: qf retraction defined in [AMS2008, (4.8)]
	POLAR: polar based retraction defined in [AMS2008, (4.7)]
	EXP: The exponential mapping
	CONSTRUCTED: the constructed retraction based on the idea in [HGA2015, Section 4.3], derivation is in [Hua2013, Section 10.2.3]
	CAYLEYR: the Cayley transform in [Zhu2016]
	[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
	Princeton University Press, Princeton, NJ, 2008.
	[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
	SIAM Journal on Optimization, 25(3):1660?685,2015.
	[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
	PhD thesis, Florida State University, Department of Mathematics, 2013.
	[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
	enum StieRetraction{ QF, POLAR, EXP, CONSTRUCTED, CAYLEYR, PROXSTIE, STIERETRACTIONLENGTH };

	/*Vector transport for the Stiefel manifold
	PARALLELIZATION: Vector transport by parallelization, See [HAG2015, Section 2.3.1]
	RIGGING: Vector transport by rigging, See [HAG2015, Section 2.3.2]
	PARALLELTRANSLATION: parallel translation
	CAYLEYVT: the vector transport based on Cayley transform. [Zhu2016]
	[HAG2015]:W. Huang, P.-A. Absil, and K. A. Gallivan. A Riemannian symmetric rank-one trust-region method.
	Mathematical Programming, 150(2):179?16, February 2015
	[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold */
	enum StieVectorTransport{ PARALLELIZATION, RIGGING, PARALLELTRANSLATION, PROJECTION, CAYLEYVT, STIEVECTORTRANSPORTLENGTH };

	class Stiefel : public Manifold{
	public:
		/*Construct the Stiefel manifold: St(p, n) and set up default parameters*/
		Stiefel(integer n, integer p);

		/*Delete EMPTYINTR and EMPTYEXTR*/
		virtual ~Stiefel(void);

		/* choose Euclidean metric, qf, parallelization and intrinsic approach and use householder reflections*/
		virtual void ChooseStieParamsSet1(void);

		/* choose Euclidean metric,  constructed retraction, parallelization and intrinsic representation */
		virtual void ChooseStieParamsSet2(void);

		/* choose Euclidean metric,  qf retraction, vector transport by projection and extrinsic representation
		TODO */
		virtual void ChooseStieParamsSet3(void);

		/* choose Euclidean metric,  Cayley retraction, Cayley vector transport and extrinsic representation
		TODO */
		virtual void ChooseStieParamsSet4(void);

		/*Euclidean metric*/
		virtual double Metric(Variable *x, Vector *etax, Vector *xix) const;

		/*Call a member function "IntrProjection" or "ExtrProjection" based on member variable "IsIntrApproach"*/
		virtual void Projection(Variable *x, Vector *v, Vector *result) const;

		/*Call a member function "qfRetraction" or "ConRetraction" based on member variable "retraction". */
		virtual void Retraction(Variable *x, Vector *etax, Variable *result) const;

		/*Call a member function "qfcoTangentVector" or "ConcoTangentVector" based on member variable "retraction". */
		virtual void coTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*Call a member function "DiffqfRetraction" or "DiffConRetraction" based on member variable "retraction". */
		virtual void DiffRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Obtain beta = \|etax\| / \|\mathcal{T}_{R_etax} etax\|
		beta has computed in "DiffRetraction". It is not necessary to recompute it in this function. */
		virtual double Beta(Variable *x, Vector *etax) const;

		/*The implementation of vector transport by parallization is identity under the intrinsic representation.
		If one needs to use the idea in [HGA2015, Section 4.1], then Manifold::LCVectorTransport is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
		virtual void VectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*The implementation of inverse vector transport by parallization is identity under the intrinsic representation.
		If one needs to use the idea in [HGA2015, Section 4.1], then Manifold::LCInverseVectorTransport is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
		virtual void InverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*The implementation of inverse vector transport by parallization is identity under the intrinsic representation.
		Therefore, Manifold::HInvTran can be used directly. If one needs to use the idea in [HGA2015, Section 4.1], then
		Manifold::LCHInvTran is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
		virtual void HInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*The implementation of vector transport by parallization is identity under the intrinsic representation.
		Therefore, Manifold::TranH can be used directly. If one needs to use the idea in [HGA2015, Section 4.1],
		then Manifold::LCTranH is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
		virtual void TranH(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, integer start, integer end, LinearOPE *result) const;

		/*The implementation of vector transport by parallization and inverse vector transport by parallelization are identity under
		the intrinsic representation. Therefore, Manifold::TranHInvTran can be used directly. If one needs to use the idea in [HGA2015, Section 4.1],
		then Manifold::LCTranHInvTran is called instead.
		[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
		SIAM Journal on Optimization, 25(3):1660?685,2015. */
		virtual void TranHInvTran(Variable *x, Vector *etax, Variable *y, LinearOPE *Hx, LinearOPE *result) const;

		/*Call a member function "ObtainIntrHHR" or "ObtainIntrSquare" based on member variable "retraction". */
		virtual void ObtainIntr(Variable *x, Vector *etax, Vector *result) const;

		/*Call a member function "ObtainExtrHHR" or "ObtainExtrSquare" based on member variable "retraction". */
		virtual void ObtainExtr(Variable *x, Vector *intretax, Vector *result) const;

		/*IntrProjection is identity, i.e., result <-- v*/
		virtual void IntrProjection(Variable *x, Vector *v, Vector *result) const;

		/*ExtrProjection is given in e.g., [Hua2013, (10.2.5)], i.e., result = v - x sym(x^T v), where sym(M) = (M + M^T)/2.
			[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
			PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		virtual void ExtrProjection(Variable *x, Vector *v, Vector *result) const;

		/*Check whether all the parameters are legal or not.*/
		virtual void CheckParams(void) const;

		/*When the metric is Euclidean, the Riemannian gradient is obtained by projecting the Euclidean
		onto the tangent space of x.*/
		virtual void EucGradToGrad(Variable *x, Vector *egf, Vector *gf, const Problem *prob) const;

		/*When the metric is Euclidean, the Riemannian action of the Hessian is obtained by
		Hess f(x)[etax] = P_x(D grad f(x) [etax]).*/
		virtual void EucHvToHv(Variable *x, Vector *etax, Vector *exix, Vector* xix, const Problem *prob) const;

		/*PARAMSMAP is defined in "def.h" and it is a map from string to double, i.e., std::map<std::string, double> .
		This function is used to set the parameters by the mapping*/
		virtual void SetParams(PARAMSMAP params);
	protected:
		integer n; /*The number of row*/
		integer p; /*The number of column*/
		StieMetric metric; /*Riemannian metric*/
		StieRetraction retraction; /*The used retraction*/
		StieVectorTransport VecTran; /*The used vector transport*/

		/*Householder transformations are used to obtain the intrinsic representation of etax*/
		virtual void ObtainIntrHHR(Variable *x, Vector *etax, Vector *result) const;

		/*Householder transformations are used to obtain the extrinsic representation of intretax*/
		virtual void ObtainExtrHHR(Variable *x, Vector *intretax, Vector *result) const;

		/*qf retraction defined in [AMS2008, (4.8)]
			[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
			Princeton University Press, Princeton, NJ, 2008.*/
		virtual void qfRetraction(Variable *x, Vector *etax, Variable *result) const;

		/*the cotangent vector for the qf retraction in [Hua2013, Section 10.2.4]
			[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
			PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		virtual void qfcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*the vector transport by differentiated retraction, see [AMS2008, Example 8.1.5]
			[AMS2008]P.-A. Absil, R. Mahony, and R. Sepulchre. Optimization algorithms on matrix manifolds.
			Princeton University Press, Princeton, NJ, 2008.*/
		virtual void DiffqfRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*the constructed retraction based on the idea in [HGA2015, Section 4.3], derivation is in [Hua2013, Section 10.2.3]
			[HGA2015]:Wen Huang, K. A. Gallivan, and P.-A. Absil. A Broyden Class of Quasi-Newton Methods for Riemannian Optimization.
			SIAM Journal on Optimization, 25(3):1660?685,2015.
			[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
			PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		virtual void ConRetraction(Variable *x, Vector *etax, Variable *result) const;

		/*the cotangent vector for the constructed retraction in [Hua2013, Section 10.2.4]
			[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
			PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		virtual void ConcoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*the vector transport by differentiated retraction, see [Hua2013, (10.2.17)]
			[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
			PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		virtual void DiffConRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*Cayley retraction defined in [Zhu2016, (15)]
		[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold.*/
		virtual void CayleyRetraction(Variable *x, Vector *etax, Variable *result) const;

		/*the cotangent vector for the Cayley retraction
		[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold.*/
		virtual void CayleycoTangentVector(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*the vector transport by differentiated retraction, see [Zhu2016, (18)]
		[Zhu2016]: Xiaojing Zhu, A Riemannian conjugate gradient method for optimization on the Stiefel Manifold.*/
		virtual void DiffCayleyRetraction(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result, bool IsEtaXiSameDir = false) const;

		/*The implementation of the isometric vector transport in [Zhu2016, (22)] */
		virtual void CayleyVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xix, Vector *result) const;

		/*The implementation of inverse vector transport, it is the inverse of [Zhu2016, (22)] */
		virtual void CayleyInverseVectorTransport(Variable *x, Vector *etax, Variable *y, Vector *xiy, Vector *result) const;

		/*compute x_perp for the constructed retraction. This is used in functions "ObtainIntrSquare" and "ObtainExtrSquare".
		The idea of computing x_perp follows from [Hua2013, Section 10.2.3]
		[Hua2013]:W. Huang. Optimization algorithms on Riemannian manifolds with applications.
		PhD thesis, Florida State University, Department of Mathematics, 2013.*/
		void ObtainPerp(Variable *x) const;

		/*x and x_perp obtained by member function "ObtainPerp" are used to obtain intrinsic representation.*/
		virtual void ObtainIntrSquare(Variable *x, Vector *etax, Vector *result) const;

		/*x and x_perp obtained by member function "ObtainPerp" are used to obtain extrinsic representation.*/
		virtual void ObtainExtrSquare(Variable *x, Vector *intretax, Vector *result) const;
	};
}; /*end of ROPTLIB namespace*/

#endif
