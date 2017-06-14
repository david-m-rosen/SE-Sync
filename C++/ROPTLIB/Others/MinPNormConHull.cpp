
#include "Others/MinPNormConHull.h"

/*Define the namespace*/
namespace ROPTLIB{

	double MinPNormConHull(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, Vector *Soln, double *YtY, integer inc)
	{
		return MinPNormConHull(Mani, x, Ys, LYs, nullptr, nullptr, Soln, YtY, inc);
	};

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y|| */
	double MinPNormConHull(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, QuasiNewton *solver, 
		void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc)
	{
#ifdef MATLAB_MEX_FILE
		return MinPNormConHullRMethod(Mani, x, Ys, LYs, solver, Hv, Soln, YtPY, inc);
//		return MinPNormConHullMatlab(Mani, x,Ys, LYs, solver, Hv, Soln, YtPY, inc);
#else
#ifdef RIEMANNIANCONHULL
		/*Use Riemannian method*/
		return MinPNormConHullRMethod(Mani, x, Ys, LYs, solver, Hv, Soln, YtPY, inc);
#else
#ifdef RECURSIVEMETHOD
		return MinPNormConHullRecursive(Mani, x, Ys, LYs, solver, Hv, Soln, YtPY, inc);
#endif 
#endif 
#endif
	};

	/*Compute min_{y in convex hull of gfs, and gfs are tangent vectors at the tangent space at x} ||y|| */
	double MinPNormConHullMatlab(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, QuasiNewton *solver, 
		void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc)
	{
		double result = 0;
#ifdef MATLAB_MEX_FILE
		integer dim = Ys[0]->Getlength();
		mxArray *GtPG = mxCreateDoubleMatrix(LYs, LYs, mxREAL);
		double *GtPGptr = mxGetPr(GtPG);

		double *G = new double[dim * LYs];
		for (integer i = 0; i < LYs; i++)
		{
			const double *tmp = Ys[i]->ObtainReadData();
			dcopy_(&dim, const_cast<double *>(tmp), &GLOBAL::IONE, G + i * dim, &GLOBAL::IONE);
		}

		if (YtPY == nullptr)
		{
			Vector **PYs = nullptr;
			if (solver != nullptr)
			{
				PYs = new Vector *[LYs];
				for (integer i = 0; i < LYs; i++)
				{
					PYs[i] = Ys[0]->ConstructEmpty();
					(solver->*Hv)(Ys[i], PYs[i]);
				}
				double *PG = new double[dim * LYs];

				for (integer j = 0; j < LYs; j++)
				{
					const double *tmp = PYs[j]->ObtainReadData();
					dcopy_(&dim, const_cast<double *>(tmp), &GLOBAL::IONE, PG + j * dim, &GLOBAL::IONE);
				}

				dgemm_(GLOBAL::T, GLOBAL::N, &LYs, &LYs, &dim, &GLOBAL::DONE, G, &dim, PG, &dim, &GLOBAL::DZERO, GtPGptr, &LYs);
				delete[] PG;
			}
			else
			{
				dgemm_(GLOBAL::T, GLOBAL::N, &LYs, &LYs, &dim, &GLOBAL::DONE, G, &dim, G, &dim, &GLOBAL::DZERO, GtPGptr, &LYs);
			}
			if (PYs != nullptr)
			{
				for (integer i = 0; i < LYs; i++)
				{
					delete PYs[i];
				}
				delete[] PYs;
			}
		}
		else
		{
			for (integer i = 0; i < LYs - 1; i++)
				for (integer j = 0; j < LYs - 1; j++)
					GtPGptr[j + i * LYs] = YtPY[j + i * inc];

			Vector *PYlast = Ys[0]->ConstructEmpty();
			(solver->*Hv)(Ys[LYs - 1], PYlast);

			dgemm_(GLOBAL::T, GLOBAL::N, &LYs, &GLOBAL::IONE, &dim, &GLOBAL::DONE, G, &dim, const_cast<double *> (PYlast->ObtainReadData()), &dim, 
				&GLOBAL::DZERO, GtPGptr + (LYs - 1) * LYs, &LYs);

			for (integer i = 0; i < LYs; i++)
			{
				GtPGptr[(LYs - 1) + i * LYs] = GtPGptr[i + (LYs - 1) * LYs];
				YtPY[(LYs - 1) + i * inc] = GtPGptr[i + (LYs - 1) * LYs];
				YtPY[i + (LYs - 1) * inc] = GtPGptr[i + (LYs - 1) * LYs];
			}

			delete PYlast;
		}

		mxArray *lhs[2], *rhs[2];
		rhs[0] = mxCreateString("ConHullProb");
		rhs[1] = GtPG;

		mexCallMATLAB(2, lhs, 2, rhs, "feval");
		if(Soln != nullptr)
		{
			double *solnptr = mxGetPr(lhs[0]);
			Mani->ScaleTimesVector(x, solnptr[0], Ys[0], Soln);
			for (integer i = 1; i < LYs; i++)
			{
				Mani->scalarVectorAddVector(x, solnptr[i], Ys[i], Soln, Soln);
			}
		}
		result = mxGetScalar(lhs[1]);

		mxDestroyArray(lhs[0]);
		mxDestroyArray(lhs[1]);
		mxDestroyArray(GtPG);
		mxDestroyArray(rhs[0]);
		delete[] G;
#endif
		return fabs(result);
	};

	double MinPNormConHullRMethod(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, QuasiNewton *solver,
		void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc)
	{
		SphereConvexHull subprob(Mani, x, Ys, LYs, solver, Hv);
		Sphere Domain(LYs);
		subprob.SetDomain(&Domain);
		SphereVariable InitX(LYs);
		InitX.RandInManifold();

		Solvers *RTRNewtonsolver = new RTRNewton(&subprob, &InitX);
		RTRNewtonsolver->Stop_Criterion = GRAD_F;
		RTRNewtonsolver->Debug = NOOUTPUT;
		RTRNewtonsolver->Max_Iteration = 100;
		RTRNewtonsolver->Tolerance = 1e-10;
		RTRNewtonsolver->Run();
		double result = RTRNewtonsolver->Getfinalfun();

		if (Soln != nullptr)
		{
			const Variable *Xopt = RTRNewtonsolver->GetXopt();
			const SharedSpace *SharedWxsq = Xopt->ObtainReadTempData("Wxsq");
			SharedWxsq->GetSharedElement()->CopyTo(Soln);
		}
		delete RTRNewtonsolver;
		return result;
	};

	double MinPNormConHullRecursive(const Manifold *Mani, Variable *x, Vector **Ys, integer LYs, QuasiNewton *solver, 
		void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, double *YtPY, integer inc)
	{
		/*The following use the method in [SY1992]*/
		/*X and Y are iterates*/
		Vector *X = Ys[0]->ConstructEmpty();
		Vector *Y = X->ConstructEmpty();
		Vector *PX = X->ConstructEmpty();
		Vector *PY = X->ConstructEmpty();
		Vector *PmY = X->ConstructEmpty();
		Vector *PYmPX = X->ConstructEmpty();

		/*The initiate for X is the one with minimum P-norm*/
		Vector **PYs = new Vector *[LYs];
		for (integer i = 0; i < LYs; i++)
		{
			PYs[i] = Ys[0]->ConstructEmpty();
		}
		if (solver != nullptr)
		{
			for (integer i = 0; i < LYs; i++)
				(solver->*Hv)(Ys[i], PYs[i]);
		}
		else
		{
			for (integer i = 0; i < LYs; i++)
				Ys[i]->CopyTo(PYs[i]);
		}

		double *Pnorm2Ys = new double[LYs];
		for (integer i = 0; i < LYs; i++)
		{
			Pnorm2Ys[i] = Mani->Metric(x, Ys[i], PYs[i]);
		}
		double minPn2 = Pnorm2Ys[0], maxPn2 = Pnorm2Ys[0];
		integer minPn2idx = 0;
		for (integer i = 1; i < LYs; i++)
		{
			if (minPn2 > Pnorm2Ys[i])
			{
				minPn2 = Pnorm2Ys[i];
				minPn2idx = i;
			}
			if (maxPn2 < Pnorm2Ys[i])
				maxPn2 = Pnorm2Ys[i];
		}
		Ys[minPn2idx]->CopyTo(X);

		double NumErr = sqrt(std::numeric_limits<double>::epsilon());//-- *LYs * ((maxPn2 > 10) ? maxPn2 : 10);

		double alpha = 0, beta = 0, lambda = 0, tmp = 0, tmp2 = 0;
		double result = 0;
		double *XtPYs = new double[LYs];
		bool *idxYs = new bool[LYs];
		integer times = 0;

		/*Start the loop*/
		while (1)
		{
			/*Step 1 in [SY1993]*/
			alpha = 1e16;
			for (integer i = 0; i < LYs; i++)
			{
				XtPYs[i] = Mani->Metric(x, X, PYs[i]);
				if (alpha > XtPYs[i])
				{
					alpha = XtPYs[i];
				}
			}
			if (solver != nullptr)
				(solver->*Hv)(X, PX);
			else
				X->CopyTo(PX);
			result = Mani->Metric(x, X, PX);
			if (result <= alpha + NumErr)
			{
				if (Soln != nullptr)
					X->CopyTo(Soln);
				break;
			}

			/*Step 2 in [SY1993]*/
			for (integer i = 0; i < LYs; i++) /*Pk is initialized by an empty set.*/
				idxYs[i] = false;

			for (integer i = 0; i < LYs; i++)
			{
				if (fabs(alpha - XtPYs[i]) <= NumErr)
				//if (alpha == XtPYs[i])
					idxYs[i] = true;
			}
			MinPNormConHullRecursivesubprob(Mani, x, Ys, PYs, LYs, idxYs, Pnorm2Ys, NumErr, Y, solver, Hv, Y, times);

			/*Step 3*/
			beta = 1e16;
			for (integer i = 0; i < LYs; i++)
			{
				if (!idxYs[i])
				{
					tmp = Mani->Metric(x, Y, PYs[i]);
					if (beta > tmp)
					{
						beta = tmp;
					}
				}
			}
			if (solver != nullptr)
				(solver->*Hv)(Y, PY);
			else
				Y->CopyTo(PY);
			result = Mani->Metric(x, Y, PY);
			if (result <= beta + NumErr)
			{
				if (Soln != nullptr)
					Y->CopyTo(Soln);
				break;
			}

			/*Step 4*/
			lambda = 1e16;
			for (integer i = 0; i < LYs; i++)
			{
				if (!idxYs[i])
				{
					Mani->VectorMinusVector(x, Ys[i], Y, PmY);
					Mani->VectorMinusVector(x, PY, PX, PYmPX);
					tmp = - Mani->Metric(x, PYmPX, PmY);
					if (tmp > 0)
					{
						tmp2 = Mani->Metric(x, PmY, PX);

						if (lambda > tmp2 / tmp && tmp2 / tmp > NumErr)
						{
							lambda = tmp2 / tmp;
						}
					}
				}
			}
			Mani->VectorLinearCombination(x, 1.0 - lambda, X, lambda, Y, X);
		}

		delete[] XtPYs;
		delete[] idxYs;
		delete[] Pnorm2Ys;
		for (integer i = 0; i < LYs; i++)
		{
			delete PYs[i];
		}
		delete[] PYs;
		delete X;
		delete Y;
		delete PX;
		delete PY;
		delete PmY;
		delete PYmPX;
		return result;
	};

	double MinPNormConHullRecursivesubprob(const Manifold *Mani, Variable *x, Vector **Ys, Vector **PYs, integer LYs, bool *idxYsFull, double *Pnorm2Ys, double NumErr, Vector *initialX,
		QuasiNewton *solver, void (QuasiNewton::*Hv)(Vector *v, Vector *result), Vector *Soln, integer &times)
	{
		Vector *X = Ys[0]->ConstructEmpty();
		Vector *Y = X->ConstructEmpty();
		Vector *PX = X->ConstructEmpty();
		Vector *PY = X->ConstructEmpty();
		Vector *PmY = X->ConstructEmpty();
		Vector *PYmPX = X->ConstructEmpty();
		times++;

		bool *idxYs = new bool[LYs];
		for (integer i = 0; i < LYs; i++)
			idxYs[i] = idxYsFull[i];

		/*Choose initial iterate X*/
		if (1)//--initialX->GetSpace() == nullptr)
		{
			double minPn2 = 1e16;
			integer minPn2idx = -1;
			for (integer i = 0; i < LYs; i++)
			{
				if (minPn2 > Pnorm2Ys[i] && idxYsFull[i])
				{
					minPn2 = Pnorm2Ys[i];
					minPn2idx = i;
				}
			}
			Ys[minPn2idx]->CopyTo(X);
		}
		else
		{
			initialX->CopyTo(X);
		}

		double *XtPYs = new double[LYs];
		double result = 0, alpha = 0, beta = 0, tmp = 0, tmp2 = 0, lambda = 0;
		/*Start the loop*/
		while (1)
		{
			/*Step 1 in [SY1993]*/
			alpha = 1e16;
			for (integer i = 0; i < LYs; i++)
			{
				if (idxYsFull[i])
				{
					XtPYs[i] = Mani->Metric(x, X, PYs[i]);
					if (alpha > XtPYs[i])
					{
						alpha = XtPYs[i];
					}
				}
			}
			if (solver != nullptr)
				(solver->*Hv)(X, PX);
			else
				X->CopyTo(PX);
			result = Mani->Metric(x, X, PX);
			if (result <= alpha + NumErr)
			{
				if (Soln != nullptr)
					X->CopyTo(Soln);
				break;
			}

			/*Step 2 in [SY1993]*/
			for (integer i = 0; i < LYs; i++) /*Pk is initialized by an empty set.*/
				idxYs[i] = false;

			for (integer i = 0; i < LYs; i++)
			{
				if (fabs(alpha - XtPYs[i]) <= NumErr && idxYsFull[i])
				//if (alpha == XtPYs[i])
					idxYs[i] = true;
			}
			MinPNormConHullRecursivesubprob(Mani, x, Ys, PYs, LYs, idxYs, Pnorm2Ys, NumErr, Y, solver, Hv, Y, times);

			/*Step 3*/
			beta = 1e16;
			for (integer i = 0; i < LYs; i++)
			{
				if (idxYsFull[i] && !idxYs[i])
				{
					tmp = Mani->Metric(x, Y, PYs[i]);
					if (beta > tmp)
					{
						beta = tmp;
					}
				}
			}
			if (solver != nullptr)
				(solver->*Hv)(Y, PY);
			else
				Y->CopyTo(PY);
			result = Mani->Metric(x, Y, PY);
			if (result <= beta + NumErr)
			{
				if (Soln != nullptr)
					Y->CopyTo(Soln);
				break;
			}

			/*Step 4*/
			lambda = 1e16;
			for (integer i = 0; i < LYs; i++)
			{
				if (idxYsFull[i] && !idxYs[i])
				{
					Mani->VectorMinusVector(x, Ys[i], Y, PmY);
					Mani->VectorMinusVector(x, PY, PX, PYmPX);
					tmp = -Mani->Metric(x, PYmPX, PmY);
					if (tmp > 0)
					{
						tmp2 = Mani->Metric(x, PmY, PX);
						if (lambda > tmp2 / tmp && tmp2 / tmp > NumErr)
						{
							lambda = tmp2 / tmp;
						}
					}
				}
			}
			Mani->VectorLinearCombination(x, 1.0 - lambda, X, lambda, Y, X);
		}
		delete[] idxYs;
		delete[] XtPYs;
		delete X;
		delete Y;
		delete PX;
		delete PY;
		delete PmY;
		delete PYmPX;
		return result;
	};

}; /*end of ROPTLIB namespace*/
