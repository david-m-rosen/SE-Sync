
#include "Problems/ElasticCurvesRO/ElasticCurvesRO.h"

/*Define the namespace*/
namespace ROPTLIB{

	ElasticCurvesRO::ElasticCurvesRO(double *inq1, double *inq2, integer ind, integer inn, double inw, bool inrotated, bool inisclosed)
	{
		n = inn;
		d = ind;
		w = inw;
		rotated = inrotated;
		isclosed = inisclosed;
		q1 = inq1;
		q2_coefs = new double[4 * d * (n - 1) + 3 * d * (n - 1) + 2 * d * (n - 1)];
		dq2_coefs = q2_coefs + 4 * d * (n - 1);
		ddq2_coefs = dq2_coefs + 3 * d * (n - 1);
		if (isclosed)
		{
			for (integer i = 0; i < d; i++)
			{
				Spline::SplineUniformPeriodic(inq2 + i * n, n, 1.0 / (n - 1), q2_coefs + i * 4 * (n - 1));
			}
		}
		else
		{
			for (integer i = 0; i < d; i++)
			{
				Spline::SplineUniformSlopes(inq2 + i * n, n, 1.0 / (n - 1), q2_coefs + i * 4 * (n - 1));
			}
		}
		for (integer i = 0; i < d; i++)
		{
			Spline::FirstDeri(q2_coefs + i * 4 * (n - 1), n, dq2_coefs + i * 3 * (n - 1));
			Spline::SecondDeri(q2_coefs + i * 4 * (n - 1), n, ddq2_coefs + i * 2 * (n - 1));
		}
	};

	ElasticCurvesRO::~ElasticCurvesRO()
	{
		delete[] q2_coefs;
	};

	double ElasticCurvesRO::f(Variable *x) const
	{
		if (x->TempDataExist("w"))
		{
			const SharedSpace *Sharedw = x->ObtainReadTempData("w");
			const double *wptr = Sharedw->ObtainReadData();
			w = wptr[0];
		}
		else
		{
			SharedSpace *Sharedw = new SharedSpace(1, 1);
			double *wptr = Sharedw->ObtainWriteEntireData();
			wptr[0] = w;
			x->AddToTempData("w", Sharedw);
		}
		const double *l = x->ObtainReadData();
		const double *O = l + n;
		const double *m = O + d * d;
		// ForDebug::Print("x:", l, 4);//---
		//double *l = x->ObtainWritePartialData();
		//double *O = l + n;
		//double *m = O + d * d;
		//O[0] = 0.877583;
		//O[1] = 0.479426;
		//O[2] = -0.479426;
		//O[3] = 0.877583;
		//x->Print("x:");//---
		//ForDebug::Print("O in f:", O, d, d);//---
		// obtain gamma
		SharedSpace *Sharedgam = new SharedSpace(1, n);
		double *gamma = Sharedgam->ObtainWriteEntireData();
		gamma[0] = 0;
		double tmp = l[0] * l[0], tmp2;
		integer deno = 2 * (n - 1);
		for (integer i = 1; i < n; i++)
		{
			tmp2 = l[i] * l[i];
			gamma[i] = gamma[i - 1] + (tmp2 + tmp) / deno;
			tmp = tmp2;
		}

		if (isclosed)
		{
			double m0 = std::fmod(m[0], 1);
			while (m0 > 1)
			{
				m0 -= 1.0;
			}
			while (m0 < 0)
			{
				m0 += 1.0;
			}
			for (integer i = 0; i < n; i++)
			{
				gamma[i] += m0;
				gamma[i] = (gamma[i] > 1) ? gamma[i] - 1 : gamma[i];
			}
		}
		// obtain q2 \circ gamma
		SharedSpace * Sharedq2g = new SharedSpace(1, n * d);
		double *q2g = Sharedq2g->ObtainWriteEntireData();
		double intv = 1.0 / (n - 1);

		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < d; j++)
			{
				q2g[i + j * n] = Spline::ValSplineUniform(q2_coefs + j * 4 * (n - 1), n, intv, gamma[i]);
			}
		}
		// obtain O q1
		SharedSpace *SharedOq1 = nullptr;
		double *Oq1 = q1;
		if (rotated)
		{
			SharedOq1 = new SharedSpace(1, n * d);
			Oq1 = SharedOq1->ObtainWriteEntireData();
			char *transn = const_cast<char *> ("n");
			char *transt = const_cast<char *> ("t");
			double one = 1, zero = 0;
			// Oq1 = q1 * O^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
			dgemm_(transn, transt, &n, &d, &d, &one, q1, &n, const_cast<double *> (O), &d, &zero, Oq1, &n);
		}

		// Oq1 - q2 l(t)
		SharedSpace *SharedOq1mq2l = new SharedSpace(1, n * d);
		double *Oq1mq2l = SharedOq1mq2l->ObtainWriteEntireData();
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < d; j++)
			{
				Oq1mq2l[i + j * n] = Oq1[i + j * n] - q2g[i + j * n] * l[i];
			}
		}
		double result = 0;
		tmp = 0;
		for (integer i = 0; i < d; i++)
		{
			tmp2 = Oq1mq2l[0 + i * n];
			tmp += tmp2 * tmp2;
		}
		result = tmp / 2;
		for (integer i = 1; i < n - 1; i++)
		{
			tmp = 0;
			for (integer j = 0; j < d; j++)
			{
				tmp2 = Oq1mq2l[i + j * n];
				tmp += tmp2 * tmp2;
			}
			result += tmp;
		}
		tmp = 0;
		for (integer j = 0; j < d; j++)
		{
			tmp2 = Oq1mq2l[n - 1 + n * j];
			tmp += tmp2 * tmp2;
		}
		result += tmp / 2;
		result /= (n - 1);

		// add penlty term
		double penlty = 0;
		tmp = l[0] * l[0];
		tmp = (tmp + 1.0 / tmp) * sqrt(1.0 + tmp * tmp);
		for (integer i = 1; i < n; i++)
		{
			tmp2 = l[i] * l[i];
			tmp2 = (tmp2 + 1.0 / tmp2) * sqrt(1.0 + tmp2 * tmp2);
			penlty += (tmp2 + tmp) / 2;
			tmp = tmp2;
		}
		penlty *= w / (n - 1);
		result += penlty;

		// attach data to x. the data can be used in gradient and hessian computation.
		x->AddToTempData("q2g", Sharedq2g);
		x->AddToTempData("gamma", Sharedgam);
		x->AddToTempData("Oq1md2l", SharedOq1mq2l);
		if (rotated)
		{
			x->AddToTempData("Oq1", SharedOq1);
		}
		return result;
	};

	void ElasticCurvesRO::EucGrad(Variable *x, Vector *egf) const
	{
		const double *l = x->ObtainReadData();

		const SharedSpace *Sharedq2g = x->ObtainReadTempData("q2g");
		const double *q2g = Sharedq2g->ObtainReadData();
		const SharedSpace *Sharedgam = x->ObtainReadTempData("gamma");
		const double *gam = Sharedgam->ObtainReadData();
		const double *Oq1 = q1;
		if (rotated)
		{
			const SharedSpace *SharedOq1 = x->ObtainReadTempData("Oq1");
			Oq1 = SharedOq1->ObtainReadData();
		}
		// obtain q2' \circ gamma 
		SharedSpace *Shareddq2g = new SharedSpace(1, n * d);
		double *dq2g = Shareddq2g->ObtainWriteEntireData();
		double intv = 1.0 / (n - 1);

		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < d; j++)
			{
				dq2g[i + j * n] = Spline::ValFirstDeriUniform(dq2_coefs + j * 3 * (n - 1), n, intv, gam[i]);
			}
		}
		// x, dy
		double *xx = new double[2 * n];
		double *dyy = xx + n;

		const SharedSpace *SharedOq1mq2l = x->ObtainReadTempData("Oq1md2l");
		const double *Oq1mq2l = SharedOq1mq2l->ObtainReadData();
		for (integer i = 0; i < n; i++)
		{
			xx[i] = 0;
			dyy[i] = 0;
			for (integer j = 0; j < d; j++)
			{
				xx[i] += Oq1mq2l[i + j * n] * q2g[i + j * n];
				dyy[i] += Oq1mq2l[i + j * n] * l[i] * dq2g[i + j * n] * 2.0;
			}
		}

		SharedSpace *Sharedyy = new SharedSpace(1, n);
		double *yy = Sharedyy->ObtainWriteEntireData();

		double *dl = egf->ObtainWriteEntireData();
		double *dO = dl + n;
		double *dm = dO + d * d;
		// compute dl
		integer deno = 2 * (n - 1);
		yy[0] = 0;
		for (integer i = 1; i < n; i++)
		{
			yy[i] = yy[i - 1] + (dyy[i] + dyy[i - 1]) / deno;
		}
		for (integer i = 0; i < n; i++)
		{
			dl[i] = 2.0 * (yy[i] * l[i] - xx[i]);
		}

		// gradeint of penlty term
		double tmp = 0, tmp2 = 0;
		for (integer i = 0; i < n; i++)
		{
			tmp = l[i] * l[i];
			tmp *= tmp;
			tmp2 = sqrt(1.0 + tmp);
			dl[i] += w * 2.0 * l[i] * (2.0 - 1.0 / tmp) * tmp2;
			//dl[i] = w * 2.0 * l[i] * (2.0 - 1.0 / tmp) * tmp2;//---
		}
		//ForDebug::Print("l:", l, n);//-----
		//ForDebug::Print("l:", l, n);//---
		//ForDebug::Print("dl:", dl, n);//---
		// compute dO;
		if (rotated)
		{
			for (integer i = 0; i < d; i++)
			{
				for (integer j = 0; j < d; j++)
				{
					dO[i + j * d] = q2g[0 + i * n] * q1[0 + j * n] * l[0] / 2;
					for (integer k = 1; k < n - 1; k++)
					{
						dO[i + j * d] += q2g[k + i * n] * q1[k + j * n] * l[k];
					}
					dO[i + j * d] += q2g[n - 1 + i * n] * q1[n - 1 + j * n] * l[n - 1] / 2;
					dO[i + j * d] *= -2.0 / (n - 1);
				}
			}
		}
		else
		{
			for (integer i = 0; i < d * d; i++)
			{
				dO[i] = 0;
			}
		}
		//ForDebug::Print("dO:", dO, d, d);//---
		// compute dm
		if (isclosed)
		{
			double tmp = 0;
			dm[0] = 0;
			for (integer j = 0; j < d; j++)
			{
				tmp += Oq1mq2l[0 + j * n] * l[0] * dq2g[0 + j * n];
			}
			tmp /= 2;
			dm[0] = tmp;
			for (integer i = 1; i < n - 1; i++)
			{
				tmp = 0;
				for (integer j = 0; j < d; j++)
				{
					tmp += Oq1mq2l[i + j * n] * l[i] * dq2g[i + j * n];
				}
				dm[0] += tmp;
			}
			tmp = 0;
			for (integer j = 0; j < d; j++)
			{
				tmp += Oq1mq2l[n - 1 + j * n] * l[n - 1] * dq2g[n - 1 + j * n];
			}
			dm[0] += tmp / 2;
			dm[0] *= -2.0 / (n - 1);
		}
		else
		{
			dm[0] = 0;
		}

		if (UseHess)
		{
			x->AddToTempData("dq2g", Shareddq2g);
			x->AddToTempData("yy", Sharedyy);
		}
		else
		{
			delete Shareddq2g;
			delete Sharedyy;
		}
		delete[] xx;
	};

	void ElasticCurvesRO::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		const double *l = x->ObtainReadData();

		const SharedSpace *Sharedq2g = x->ObtainReadTempData("q2g");
		const double *q2g = Sharedq2g->ObtainReadData();
		const SharedSpace *Sharedgam = x->ObtainReadTempData("gamma");
		const double *gam = Sharedgam->ObtainReadData();
		const double *Oq1 = q1;
		if (rotated)
		{
			const SharedSpace *SharedOq1 = x->ObtainReadTempData("Oq1");
			Oq1 = SharedOq1->ObtainReadData();
		}
		const SharedSpace *SharedOq1mq2l = x->ObtainReadTempData("Oq1md2l");
		const double *Oq1mq2l = SharedOq1mq2l->ObtainReadData();
		const SharedSpace *Shareddq2g = x->ObtainReadTempData("dq2g");
		const double *dq2g = Shareddq2g->ObtainReadData();
		const SharedSpace *Sharedyy = x->ObtainReadTempData("yy");
		const double *yy = Sharedyy->ObtainReadData();

		const double *v = etax->ObtainReadData();
		const double *U = v + n;
		const double *a = U + d * d;

		double *Uq1 = new double[n * d + 2 * n + 3 * n * d + 2 * n];
		double *tmpn = Uq1 + d * n;
		double *tmpn2 = tmpn + n;
		double *Dq2gva = tmpn2 + n;
		double *Ddq2gva = Dq2gva + d * n;
		double *tmpnd = Ddq2gva + d * n;
		double *Dx = tmpnd + d * n;
		double *Dy = Dx + n; // n

		// compute Dq2gva and Ddq2gva
		PointwiseProd(v, l, n, tmpn);
		CumTrapz(tmpn, n, 1.0 / (n - 1), tmpn2);
		PointwiseQProdl(dq2g, tmpn2, d, n, Dq2gva);
		integer nd = n * d, inc = 1;
		double one = 1;
		// Dq2gva <- a * dq2g + Dq2gva, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&nd, const_cast<double *> (a), const_cast<double *> (dq2g), &inc, Dq2gva, &inc);

		if (!x->TempDataExist("ddq2g"))
		{
			SharedSpace *Sharedddq2g = new SharedSpace(1, n * d);
			double *ddq2g = Sharedddq2g->ObtainWriteEntireData();
			for (integer i = 0; i < n; i++)
			{
				for (integer j = 0; j < d; j++)
				{
					ddq2g[i + j * n] = Spline::ValSecondDeriUniform(ddq2_coefs + j * 2 * (n - 1), n, 1.0 / (n - 1), gam[i]);
				}
			}
			x->AddToTempData("ddq2g", Sharedddq2g);
		}
		const SharedSpace *Sharedddq2g = x->ObtainReadTempData("ddq2g");
		const double *ddq2g = Sharedddq2g->ObtainReadData();

		PointwiseQProdl(ddq2g, tmpn2, d, n, Ddq2gva);
		// Ddq2gva <- a * ddq2g + Ddq2gva, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&nd, const_cast<double *> (a), const_cast<double *> (ddq2g), &inc, Ddq2gva, &inc);

		// compute DX
		char *transn = const_cast<char *> ("n");
		double zero = 0;
		// Uq1 = q1 * U, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &n, &d, &d, &one, q1, &n, const_cast<double *> (U), &d, &zero, Uq1, &n);
		PointwiseInnerProd(Uq1, q2g, d, n, tmpn);
		PointwiseInnerProd(q2g, q2g, d, n, tmpn2);
		PointwiseProd(tmpn2, v, n, tmpn2);
		double coef = -1.0;
		// tmpn <- coef * tmpn2 + tmpn, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &coef, tmpn2, &inc, tmpn, &inc);

		PointwiseQProdl(q2g, l, d, n, tmpnd);
		coef = -2.0;
		// tmpnd <- coef * tmpnd, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
		dscal_(&nd, &coef, tmpnd, &inc);
		// tmpnd <- Oq1 + tmpnd, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&nd, &one, const_cast<double *> (Oq1), &inc, tmpnd, &inc);

		PointwiseInnerProd(tmpnd, Dq2gva, d, n, Dx);
		// Dx <- tmpn + Dx, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &one, tmpn, &inc, Dx, &inc);

		// compute DY
		PointwiseInnerProd(Uq1, dq2g, d, n, Dy);
		PointwiseProd(Dy, l, n, Dy);
		double two = 2;
		// Dy <- 2 * Dy, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
		dscal_(&n, &two, Dy, &inc);

		PointwiseInnerProd(q2g, dq2g, d, n, tmpn);
		PointwiseProd(tmpn, l, n, tmpn);
		PointwiseProd(tmpn, v, n, tmpn);
		// tmpn <- 2 * tmpn, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
		dscal_(&n, &two, tmpn, &inc);
		coef = -1.0;
		// Dy <- coef * tmpn + Dy, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &coef, tmpn, &inc, Dy, &inc);

		PointwiseInnerProd(Oq1mq2l, dq2g, d, n, tmpn);
		PointwiseProd(tmpn, v, n, tmpn);
		// Dy <- 2 * tmpn + Dy, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &two, tmpn, &inc, Dy, &inc);

		PointwiseInnerProd(Dq2gva, dq2g, d, n, tmpn);
		PointwiseProd(tmpn, l, n, tmpn);
		PointwiseProd(tmpn, l, n, tmpn);
		coef = -2.0;
		// Dy <- coef * tmpn + Dy, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &coef, tmpn, &inc, Dy, &inc);

		PointwiseInnerProd(Ddq2gva, Oq1mq2l, d, n, tmpn);
		PointwiseProd(tmpn, l, n, tmpn);
		coef = 2.0;
		// Dy <- coef * tmpn + Dy, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &coef, tmpn, &inc, Dy, &inc);
		// tmpn <- Dy, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
		dcopy_(&n, Dy, &inc, tmpn, &inc);
		CumTrapz(tmpn, n, 2.0 / (n - 1), Dy); // 1.0 / (n - 1)

		double *vexix = exix->ObtainWriteEntireData();
		double *Uexix = vexix + n;
		double *mexix = Uexix + d * d;

		// compute vexix
		PointwiseProd(v, yy, n, vexix);
		PointwiseProd(l, Dy, n, tmpn);
		// vexix <- tmpn + vexix, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &one, tmpn, &inc, vexix, &inc);
		coef = -1.0;
		// vexix <- coef * Dx + vexix, details: http://www.netlib.org/lapack/explore-html/d9/dcd/daxpy_8f.html
		daxpy_(&n, &coef, Dx, &inc, vexix, &inc);
		coef = 2.0;
		// vexix <- coef * vexix, details: http://www.netlib.org/lapack/explore-html/d4/dd0/dscal_8f.html
		dscal_(&n, &coef, vexix, &inc);

		// action of Hessian of penlty term 
		double tmp = 0, tmp2 = 0;
		for (integer i = 0; i < n; i++)
		{
			tmp = l[i] * l[i];
			tmp = tmp * tmp;
			tmp2 = sqrt(1.0 + tmp);
			vexix[i] += w * 2.0 * ((2.0 + 3.0 / tmp) * tmp2 + 2.0 * (2.0 * tmp - 1.0) / tmp2) * v[i];
			//vexix[i] = w * 2.0 * ((2.0 + 3.0 / tmp) * tmp2 + 2.0 * (2.0 * tmp - 1.0) / tmp2) * v[i];//---
		}

		// compute Uexix
		if (rotated)
		{
			for (integer i = 0; i < d; i++)
			{
				for (integer j = 0; j < d; j++)
				{
					PointwiseProd(Dq2gva + i * n, q1 + j * n, n, tmpn);
					PointwiseProd(tmpn, l, n, tmpn);
					Uexix[i + j * d] = Trapz(tmpn, n, 1.0 / (n - 1));
					PointwiseProd(q2g + i * n, q1 + j * n, n, tmpn);
					PointwiseProd(tmpn, v, n, tmpn);
					Uexix[i + j * d] += Trapz(tmpn, n, 1.0 / (n - 1));
					Uexix[i + j * d] *= -2.0;
				}
			}
		}
		else
		{
			for (integer i = 0; i < d * d; i++)
				Uexix[i] = 0;
		}

		// compute mexix
		if (isclosed)
		{
			PointwiseInnerProd(Uq1, dq2g, d, n, tmpn);
			PointwiseProd(tmpn, l, n, tmpn);
			mexix[0] = Trapz(tmpn, n, 1.0 / (n - 1));

			PointwiseInnerProd(q2g, dq2g, d, n, tmpn);
			PointwiseProd(tmpn, l, n, tmpn);
			PointwiseProd(tmpn, v, n, tmpn);
			mexix[0] -= Trapz(tmpn, n, 1.0 / (n - 1));

			PointwiseInnerProd(Dq2gva, dq2g, d, n, tmpn);
			PointwiseProd(tmpn, l, n, tmpn);
			PointwiseProd(tmpn, l, n, tmpn);
			mexix[0] -= Trapz(tmpn, n, 1.0 / (n - 1));

			PointwiseInnerProd(Oq1mq2l, dq2g, d, n, tmpn);
			PointwiseProd(tmpn, v, n, tmpn);
			mexix[0] += Trapz(tmpn, n, 1.0 / (n - 1));

			PointwiseInnerProd(Oq1mq2l, Ddq2gva, d, n, tmpn);
			PointwiseProd(tmpn, l, n, tmpn);
			mexix[0] += Trapz(tmpn, n, 1.0 / (n - 1));

			mexix[0] *= -2.0;
		}
		else
		{
			mexix[0] = 0;
		}

		//etax->CopyTo(exix);
		delete[] Uq1;
	};

	void ElasticCurvesRO::PointwiseInnerProd(const double *q1, const double *q2, integer d, integer n, double *result)
	{
		for (integer i = 0; i < n; i++)
		{
			result[i] = 0;
			for (integer j = 0; j < d; j++)
			{
				result[i] += q1[i + j * n] * q2[i + j * n];
			}
		}
	};

	void ElasticCurvesRO::PointwiseQProdl(const double *q1, const double *l, integer d, integer n, double *result)
	{
		for (integer i = 0; i < n; i++)
		{
			for (integer j = 0; j < d; j++)
			{
				result[i + j * n] = q1[i + j * n] * l[i];
			}
		}
	};

	void ElasticCurvesRO::PointwiseProd(const double *l1, const double *l2, integer n, double *result)
	{
		for (integer i = 0; i < n; i++)
			result[i] = l1[i] * l2[i];
	};

	void ElasticCurvesRO::CumTrapz(const double *l, integer n, double intv, double *result)
	{
		result[0] = 0;
		for (integer i = 1; i < n; i++)
		{
			result[i] = result[i - 1] + (l[i - 1] + l[i]) * intv / 2;
		}
	};

	double ElasticCurvesRO::Trapz(const double *l, integer n, double intv)
	{
		double result = l[0] / 2;
		for (integer i = 1; i < n - 1; i++)
		{
			result += l[i];
		}
		result += l[n - 1] / 2;
		result *= intv;
		return result;
	};
}; /*end of ROPTLIB namespace*/
