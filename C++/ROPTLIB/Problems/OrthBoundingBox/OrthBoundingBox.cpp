
#include "Problems/OrthBoundingBox/OrthBoundingBox.h"

/*Define the namespace*/
namespace ROPTLIB{

	OrthBoundingBox::OrthBoundingBox(double *inE, integer ind, integer inn)
	{
		E = inE;
		d = ind;
		n = inn;
	};

	OrthBoundingBox::~OrthBoundingBox(void)
	{
	};

	double OrthBoundingBox::f(Variable *x) const
	{
		const double *xptr = x->ObtainReadData();
		Vector *gf = Domain->GetEMPTYEXTR()->ConstructEmpty();
		SharedSpace *Sharedgf = new SharedSpace(gf);
		double *gfptr = gf->ObtainWriteEntireData();
		double result = 1;

		double *XE = new double[d * n];
		Matrix MX(xptr, d, d), ME(E, d, n), MXE(XE, d, n);
		Matrix::DGEMM(1, MX, false, ME, false, 0, MXE);

		double *minv = new double[2 * d];
		double *maxv = minv + d;
		integer *minidx = new integer[2 * d];
		integer *maxidx = minidx + d;
		for (integer i = 0; i < d; i++)
		{
			minv[i] = XE[i];
			minidx[i] = 0;
			maxv[i] = XE[i];
			maxidx[i] = 0;
			for (integer j = 0; j < n; j++)
			{
				if (XE[i + j * d] < minv[i])
				{
					minv[i] = XE[i + j * d];
					minidx[i] = j;
				}
				if (XE[i + j * d] > maxv[i])
				{
					maxv[i] = XE[i + j * d];
					maxidx[i] = j;
				}
			}
			result *= (maxv[i] - minv[i]);
		}


		double tmp = 0;
		for (integer i = 0; i < d; i++)
		{
			tmp = result / (maxv[i] - minv[i]);
			for (integer j = 0; j < d; j++)
			{
				gfptr[i + j * d] = tmp * (E[j + maxidx[i] * d] - E[j + minidx[i] * d]);
			}
		}

		delete[] minidx;
		delete[] minv;
		delete[] XE;
		x->AddToTempData("gf", Sharedgf);
		return result;
	};

	void OrthBoundingBox::EucGrad(Variable *x, Vector *egf) const
	{
		const SharedSpace *Sharedgf = x->ObtainReadTempData("gf");
		Vector *gf = Sharedgf->GetSharedElement();
		gf->CopyTo(egf);
	};

	void OrthBoundingBox::AllDirectionDerivative(const Variable *x, Vector **dgfs, integer UpperBound, double threshold, integer &idxgf)
	{
		integer *maxidx = new integer[d * 2];
		integer *minidx = maxidx + d;
		double *maxv = new double[d * 2];
		double *minv = maxv + d;
		integer idx = 0;
		bool findmax = true;
		double volumn = 1;

		const double *xptr = x->ObtainReadData();
		double *XE = new double[d * n];
		Matrix MX(xptr, d, d), ME(E, d, n), MXE(XE, d, n);
		Matrix::DGEMM(1, MX, false, ME, false, 0, MXE);

		//ForDebug::Print("XE:", XE, d, n);//---

		RecursiveDirDeri(x, maxidx, minidx, maxv, minv, idx, volumn, findmax, idxgf, dgfs, UpperBound, XE, threshold);

		delete[] XE;
		delete[] maxidx;
		delete[] maxv;
	};

	void OrthBoundingBox::RecursiveDirDeri(const Variable *x, integer *maxidx, integer *minidx, double *maxv, double *minv,
		integer idx, double volumn, bool findmax, integer &idxgf, Vector **dgfs, integer UpperBound, double *XE, double threshold)
	{
		if (findmax && idx == d)
		{
			if (idxgf >= UpperBound)
			{
				idxgf++;
				return;
			}

			Vector *egf = Domain->GetEMPTYEXTR()->ConstructEmpty();
			double *gfptr = egf->ObtainWriteEntireData();
			double tmp = 0;
			for (integer i = 0; i < d; i++)
			{
				tmp = volumn / (maxv[i] - minv[i]);
				for (integer j = 0; j < d; j++)
				{
					gfptr[i + j * d] = tmp * (E[j + maxidx[i] * d] - E[j + minidx[i] * d]);
				}
			}
			Domain->ObtainIntr(const_cast<Variable *> (x), egf, dgfs[idxgf]);
			idxgf++;
			delete egf;
			return;
		}

		if (findmax)
		{
			integer imaxi = 0;
			double imaxv = XE[idx];
			/*find max value and its index.*/
			for (integer j = 0; j < n; j++)
			{
				if (XE[idx + j * d] > imaxv)
				{
					imaxv = XE[idx + j * d];
					imaxi = j;
				}
			}
			/*if the difference are with threshold, then they are considered to be max values*/
			for (integer j = 0; j < n; j++)
			{
				if (fabs(XE[idx + j * d] - imaxv) <= threshold)
				{
					maxv[idx] = XE[idx + j * d];
					maxidx[idx] = j;
					RecursiveDirDeri(x, maxidx, minidx, maxv, minv, idx, volumn, false, idxgf, dgfs, UpperBound, XE, threshold);
				}
			}
		}

		if (!findmax)
		{
			integer imini = 0;
			double iminv = XE[idx];
			/*find min value and its index.*/
			for (integer j = 0; j < n; j++)
			{
				if (XE[idx + j * d] < iminv)
				{
					iminv = XE[idx + j * d];
					imini = j;
				}
			}
			/*if the difference are with threshold, then they are considered to be min values*/
			for (integer j = 0; j < n; j++)
			{
				if (fabs(XE[idx + j * d] - iminv) <= threshold)
				{
					minv[idx] = XE[idx + j * d];
					minidx[idx] = j;
					RecursiveDirDeri(x, maxidx, minidx, maxv, minv, idx + 1, volumn * (maxv[idx] - minv[idx]), true, idxgf, dgfs, UpperBound, XE, threshold);
				}
			}
		}
	};

	void OrthBoundingBox::EucHessianEta(Variable *x, Vector *etax, Vector *exix) const
	{
		etax->CopyTo(exix);
	};
}; /*end of ROPTLIB namespace*/
