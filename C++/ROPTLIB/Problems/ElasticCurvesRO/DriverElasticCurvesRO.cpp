
#include "Problems/ElasticCurvesRO/DriverElasticCurvesRO.h"

/*Define the namespace*/
namespace ROPTLIB{

	void DriverElasticCurvesRO(double *C1, double *C2, integer d, integer n, double w, bool rotated, bool isclosed,
		bool onlyDP, integer skipm, std::string solverstr, integer autoselectC, ProductElement *Xopt, bool &swap, double *fopts, 
		double *comtime, integer &Nsout, integer &numinitialx, double *optQ1, double *optQ2)
	{ // The first and last point of C1 and C2 should be the same if they are viewed as closed curves, i.e., isclosed = true.
		double threshold = PI / 2;
		integer minSkip = skipm;
		integer randshift = 0;
		bool computeCD1 = false;

		Solvers *solver = nullptr;
		// Let C2 be the complex one
		double TAC1 = ComputeTotalAngle(C1, d, n);
		double TAC2 = ComputeTotalAngle(C2, d, n);
		double *temppt, TACtemp;
		swap = false;

		// autoselectC: 0: keep the order, 1 and 2: C2 is the simple one
		// C2 being a simple one would be prefered. The spline on a smoother may be more accurate
		// and resample on a smoother curve yields a smaller difference.
		// Points on the more complex C1 would preserve C1's structure.
		if (autoselectC != 0)
		{
			if (TAC1 < TAC2)
			{
				temppt = C1;
				C1 = C2;
				C2 = temppt;
				temppt = optQ1;
				optQ1 = optQ2;
				optQ2 = temppt;
				TACtemp = TAC1;
				TAC1 = TAC2;
				TAC2 = TACtemp;
				swap = true;
			}
		}

		// find initial breaks and Ns
		integer *ms = new integer[n];
		integer lms = 0, ns = n;
		if (isclosed)
		{
			if (onlyDP)
			{
				skipm = (skipm < 1) ? 1 : skipm;
				for (integer i = 0; i < n - 1; i += skipm)
				{
					ms[lms] = i;
					lms++;
				}
			}
			else
			{
				if (autoselectC != 1)
				{ // use complex one to generate break points, this yields a better distance estimation.
					if (TAC2 > TAC1)
					{
						FindInitialBreaksAndNs(C2, d, n, minSkip, threshold, randshift, ms, lms, ns);
						ns = static_cast<int> (static_cast<double> (n) / 3);
						ns = (ns > 30) ? 30 : ns;
						ns += static_cast<int> (TAC1 / PI * 2.0);
					}
					else
					{
						FindInitialBreaksAndNs(C1, d, n, minSkip, threshold, randshift, ms, lms, ns);
						for (integer i = 1; i < lms; i++)
						{
							ms[i] = n - ms[i];
						}
						ns = static_cast<int> (static_cast<double> (n) / 3);
						ns = (ns > 30) ? 30 : ns;
						ns += static_cast<int> (TAC1 / PI * 2.0);
					}
				}
				else
				{ // use simple one to generate break points, this yields a faster algorithm.
					if (TAC2 < TAC1)
					{
						FindInitialBreaksAndNs(C2, d, n, minSkip, threshold, randshift, ms, lms, ns);
						ns = static_cast<int> (static_cast<double> (n) / 3);
						ns = (ns > 30) ? 30 : ns;
						ns += static_cast<int> (TAC1 / PI * 2.0);
					}
					else
					{
						FindInitialBreaksAndNs(C1, d, n, minSkip, threshold, randshift, ms, lms, ns);
						for (integer i = 1; i < lms; i++)
						{
							ms[i] = n - ms[i];
						}
						ns = static_cast<int> (static_cast<double> (n) / 3);
						ns = (ns > 30) ? 30 : ns;
						ns += static_cast<int> (TAC2 / PI * 2.0);
					}
				}
			}
		}
		else
		{ // for open curve, only one break point is used.
			ms[0] = 0;
			lms = 1;
			if (!onlyDP)
			{
				ns = static_cast<int> (static_cast<double> (n) / 3);
				ns = (ns > 30) ? 30 : ns;
				ns += static_cast<int> (TAC2 / PI * 2.0);
			}
		}

		Nsout = ns;
		numinitialx = lms;

		// create manifold and initial iterate objects.
		integer numofmanis = 3;
		integer numofmani1 = 1;
		integer numofmani2 = 1;
		integer numofmani3 = 1;
		L2Sphere TNS(n);
		OrthGroup OG(d);
		Euclidean Euc(1);
		ProductManifold *Domain = nullptr;
		Domain = new ProductManifold(numofmanis, &TNS, numofmani1, &OG, numofmani2, &Euc, numofmani3);

		// 	Domain->SetIsIntrApproach(false);

		L2SphereVariable TNSV(n);
		OrthGroupVariable OGV(d);
		EucVariable EucV(1);
		ProductElement *InitialX = nullptr;
		InitialX = new ProductElement(numofmanis, &TNSV, numofmani1, &OGV, numofmani2, &EucV, numofmani3);
		double *Xptr = InitialX->ObtainWriteEntireData();

		// initialize rotation to be identity and shift to be zero
		Xptr[n + d * d] = 0;
		for (integer j = 0; j < d; j++)
		{
			Xptr[n + j + j * d] = 1;
			for (integer k = j + 1; k < d; k++)
			{
				Xptr[n + k + j * d] = 0;
				Xptr[n + j + k * d] = 0;
			}
		}


		// find initialX for each break and run the solver
		ElasticCurvesRO *ECRO = nullptr;
		double *C2shift = new double[5 * d * n + d * d + n + lms + 2 * d * d];
		double *q2shift = C2shift + d * n;
		double *q1 = q2shift + d * n;
		double *O = q1 + d * n;
		double *Rotq2shift = O + d * d;
		double *RotC2shift = Rotq2shift + d * n;
		double *DPgam = RotC2shift + d * n;
		double *msV = DPgam + n;
		double *O2 = msV + lms;
		double *O3 = O2 + d * d; // d * d

		double *C1s = nullptr, *C2s = nullptr, *q1s = nullptr, *q2s = nullptr, *DPgams = nullptr;
		if (!onlyDP)
		{
			C1s = new double[4 * d * ns + ns];
			C2s = C1s + d * ns;
			q1s = C2s + d * ns;
			q2s = q1s + d * ns;
			DPgams = q2s + d * ns; // ns
		}
		double *C2_coefs = nullptr, *q2 = nullptr;
		if (onlyDP)
		{
			C2_coefs = new double[4 * d * (n - 1) + n * d];
			q2 = C2_coefs + 4 * d * (n - 1);
		}

		if (optQ1 == nullptr)
			CurveToQ(C1, d, n, q1, isclosed);
		else
		{
			integer len = d * n;
			dcopy_(&len, optQ1, &GLOBAL::IONE, q1, &GLOBAL::IONE);
		}

		char *transn = const_cast<char *> ("n"), *transt = const_cast<char *> ("t");
		double one = 1, zero = 0;
		integer dd = d * d, inc = 1;
		unsigned long starttime = getTickCount();
		double minmsV = 10000;

		if (!onlyDP)
		{ // use coarse point to represent C1
			GetCurveSmall(C1, C1s, d, n, ns, isclosed);
			CurveToQ(C1s, d, ns, q1s, isclosed);
		}
		printf("lms:%d, ns:%d\n", lms, ns);

		double *Xoptptr = Xopt->ObtainWriteEntireData();
		Xoptptr[n + d * d] = 0;

		for (integer i = 0; i < 5; i++)
		{
			fopts[i] = 1000;
			comtime[i] = static_cast<double>(getTickCount() - starttime) / CLK_PS;
		}

		for (integer i = 0; i < lms; i++) //lms
		{
			printf("%d, ", ms[i]);
			starttime = getTickCount();
			// obtain initial reparameterization
			ShiftC(C2, d, n, C2shift, ms[i]);
			if (optQ2 == nullptr)
				CurveToQ(C2shift, d, n, q2shift, isclosed);
			else
				ShiftC(optQ2, d, n, q2shift, ms[i]);

			if (rotated)
			{ // if rotation is considered, rotate C2.
				FindBestRotation(q1, q2shift, d, n, O);
				// Rotq2shift <- q2shift * O^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
				dgemm_(transn, transt, &n, &d, &d, &one, q2shift, &n, O, &d, &zero, Rotq2shift, &n);
				// RotC2shift <- C2shift * O^T, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
				dgemm_(transn, transt, &n, &d, &d, &one, C2shift, &n, O, &d, &zero, RotC2shift, &n);
			}
			else
			{ // Otherwise, keep C2.
				integer nd = n * d, inc = 1;
				// Rotq2shift <- q2shift, details:http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
				dcopy_(&nd, q2shift, &inc, Rotq2shift, &inc);
				// RotC2shift <- C2shift, details:http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
				dcopy_(&nd, C2shift, &inc, RotC2shift, &inc);
			}

			if (!onlyDP)
			{ // if use Riemannian method, use coarse points to represent C2 and compute the initial gamma.
				GetCurveSmall(RotC2shift, C2s, d, n, ns, isclosed);
				CurveToQ(C2s, d, ns, q2s, isclosed);
				DynamicProgramming(q2s, q1s, d, ns, DPgams, isclosed, NUMSMALL);
				ReSampleGamma(DPgams, ns, DPgam, n);
			}
			else
			{ // if only use DP, then use dense grid to compute gamma
				DynamicProgramming(Rotq2shift, q1, d, n, DPgam, isclosed, NUMBIG);

				if (rotated)
				{ // if rotation is consider, then an extra rotation is computed.
					if (computeCD1)
					{ // if CD1 function value needs to be output, then compute it and output it.
						if (isclosed)
							GradientPeriod(DPgam, n, 1.0 / (n - 1), Xptr);
						else
							Gradient(DPgam, n, 1.0 / (n - 1), Xptr);
						for (integer j = 0; j < n; j++)
						{
							Xptr[j] = sqrt(Xptr[j]);
						}
						ECRO = new ElasticCurvesRO(q1, Rotq2shift, d, n, w, rotated, isclosed);
						ECRO->SetDomain(Domain);
						printf("CD1 func:%g\n", ECRO->f(InitialX));
					}

					if (isclosed)
					{
						for (integer j = 0; j < d; j++)
						{
							Spline::SplineUniformPeriodic(RotC2shift + j * n, n, 1.0 / (n - 1), C2_coefs + j * 4 * (n - 1));
						}
					}
					else
					{
						for (integer j = 0; j < d; j++)
						{
							Spline::SplineUniformSlopes(RotC2shift + j * n, n, 1.0 / (n - 1), C2_coefs + j * 4 * (n - 1));
						}
					}
					for (integer j = 0; j < n; j++)
					{
						for (integer k = 0; k < d; k++)
						{
							RotC2shift[j + k * n] = Spline::ValSplineUniform(C2_coefs + k * 4 * (n - 1), n, 1.0 / (n - 1), DPgam[j]);
						}
					}
					CurveToQ(RotC2shift, d, n, q2, isclosed);
					FindBestRotation(q1, q2, d, n, O2);
					// O3 = O * O2, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
					dgemm_(transn, transn, &d, &d, &d, &one, O, &d, O2, &d, &zero, O3, &d);
					// O <- O3, details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
					dcopy_(&dd, O3, &inc, O, &inc);
					// Xptr(n : n + d * d - 1) <- reshape(O2, d * d, 1), 
					// details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
					dcopy_(&dd, O2, &inc, Xptr + n, &inc); // used to evaluate the cost function
				}
			}

			if (isclosed)
				GradientPeriod(DPgam, n, 1.0 / (n - 1), Xptr);
			else
				Gradient(DPgam, n, 1.0 / (n - 1), Xptr);

			for (integer j = 0; j < n; j++)
			{
				Xptr[j] = sqrt(Xptr[j]);
			}
			// Compute reparameterization for q1 and rotated and shifted q2;
			ECRO = new ElasticCurvesRO(q1, Rotq2shift, d, n, w, rotated, isclosed);
			ECRO->SetDomain(Domain);
			if (onlyDP)
			{ // if only DP is used, then output the CD1H cost function
				ECRO->w = 0;
				msV[i] = ECRO->f(InitialX);
				printf("CD1H func:%g\n", msV[i]);
			}
			if (!onlyDP)
			{ // if a Riemannian method is used, then Xinitial is the initial iterate and a method is used.
				if (solverstr == "RBFGS")
				{
					solver = new RBFGS(ECRO, InitialX);
					dynamic_cast<SolversLS *> (solver)->Initstepsize = 0.001;
				}
				else
				if (solverstr == "LRBFGS")
				{
					solver = new LRBFGS(ECRO, InitialX);
					dynamic_cast<SolversLS *> (solver)->Initstepsize = 0.001;
				}
				else
				if (solverstr == "RCG")
				{
					solver = new RCG(ECRO, InitialX);
					dynamic_cast<SolversLS *> (solver)->Initstepsize = 0.001;
				}
				else
				if (solverstr == "RSD")
				{
					solver = new RSD(ECRO, InitialX);
					dynamic_cast<SolversLS *> (solver)->Initstepsize = 0.001;
				}
				else
				if (solverstr == "RTRSR1")
				{
					solver = new RTRSR1(ECRO, InitialX);
					dynamic_cast<SolversTR *> (solver)->kappa = 0.1;
					dynamic_cast<SolversTR *> (solver)->theta = 1.0;
				}
				else
				if (solverstr == "LRTRSR1")
				{
					solver = new LRTRSR1(ECRO, InitialX);
					dynamic_cast<SolversTR *> (solver)->kappa = 0.1;
					dynamic_cast<SolversTR *> (solver)->theta = 1.0;
				}
				else
				if (solverstr == "RTRSD")
				{
					solver = new RTRSD(ECRO, InitialX);
				}
				else
				{
					printf("This solver is not used in this problem!\n");
					delete ECRO;
					delete solver;
					delete[] C2shift;
					if (C2_coefs != nullptr)
					{
						delete[] C2_coefs;
					}
					if (C1s != nullptr)
					{
						delete[] C1s;
					}
					delete[] ms;
					delete Domain;
					delete InitialX;
					return;
				}
				solver->Max_Iteration = 500;
				solver->Min_Iteration = 10;
				solver->Debug = NOOUTPUT; //--FINALRESULT;//--NOOUTPUT; //ITERRESULT
				solver->Stop_Criterion = FUN_REL;
				solver->Tolerance = 1e-3;
				solver->Run();
				ECRO->w = 0;
				//--Xopt->RemoveAllFromTempData();
				msV[i] = ECRO->f(const_cast<Element *> (solver->GetXopt()));
				printf("%s func:%g, num of iter:%d\n", solverstr.c_str(), msV[i], solver->GetIter());
			}
			delete ECRO;

			if (msV[i] < minmsV)
			{
				minmsV = msV[i];
				if (onlyDP)
				{
					for (integer j = 0; j < n; j++)
					{
						Xoptptr[j] = DPgam[j];
					}
					// Xoptptr(n:n+d*d - 1) <- reshape(O, d * d, 1) 
					// details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
					dcopy_(&dd, O, &inc, Xoptptr + n, &inc);
					Xoptptr[n + d * d] = static_cast<double> (ms[i]) / (n - 1);
				}
				else
				{
					solver->GetXopt()->CopyTo(Xopt);
					Xoptptr = Xopt->ObtainWritePartialData();
					for (integer j = 0; j < n; j++)
					{
						Xoptptr[j] *= Xoptptr[j];
					}
					double tmp1 = Xoptptr[0], tmp2 = 0;
					Xoptptr[0] = 0;
					for (integer j = 1; j < n; j++)
					{
						tmp2 = Xoptptr[j];
						Xoptptr[j] = Xoptptr[j - 1] + (tmp1 + tmp2) / 2 / (n - 1);
						tmp1 = tmp2;
					}
					// O2 = O * reshape(Xoptptr(n : n + d * d - 1), d, d)^T, 
					// details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
					dgemm_(transn, transt, &d, &d, &d, &one, O, &d, Xoptptr + n, &d, &zero, O2, &d);
					// Xoptptr(n:n + d * d - 1) <- reshape(O2, d * d, 1) 
					// details: http://www.netlib.org/lapack/explore-html/da/d6c/dcopy_8f.html
					dcopy_(&dd, O2, &inc, Xoptptr + n, &inc);

					Xoptptr[n + d * d] = Xoptptr[n + d * d] + static_cast<double> (ms[i]) / (n - 1);
				}
			}
			if (!onlyDP)
				delete solver;

			comtime[0] += (double)(getTickCount() - starttime) / CLK_PS;
			if (msV[i] < fopts[0])
				fopts[0] = msV[i];
			if (i % 2 == 0)
			{
				comtime[1] += (double)(getTickCount() - starttime) / CLK_PS;
				if (msV[i] < fopts[1])
					fopts[1] = msV[i];
			}
			if (i % 4 == 0)
			{
				comtime[2] += (double)(getTickCount() - starttime) / CLK_PS;
				if (msV[i] < fopts[2])
					fopts[2] = msV[i];
			}
			if (i % 8 == 0)
			{
				comtime[3] += (double)(getTickCount() - starttime) / CLK_PS;
				if (msV[i] < fopts[3])
					fopts[3] = msV[i];
			}
			if (i % 16 == 0)
			{
				comtime[4] += (double)(getTickCount() - starttime) / CLK_PS;
				if (msV[i] < fopts[4])
					fopts[4] = msV[i];
			}
		}

		printf("min f:%3.2e\n", minmsV);
		printf("time:%3.2e\n", comtime[0]);
		delete[] C2shift;
		if (C2_coefs != nullptr)
		{
			delete[] C2_coefs;
		}
		if (C1s != nullptr)
		{
			delete[] C1s;
		}
		delete[] ms;
		delete InitialX;
		delete Domain;
	};

	double DynamicProgramming(const double *q1, const double *q2, integer d, integer n, double *gamma, bool isclosed, SLOPESTYPE Nbrstype)
	{
		integer k, l, m, Eidx, Fidx, Ftmp, Fmin, Num, *Path, *x, *y, cnt;
		double *q2L, *E, Etmp, Emin, a, b;
		integer nn = n - 1, splinestatus;
		m = 5 * nn + 1;
		integer mm = m - 1;
		double mdn = static_cast<double> (mm) / (nn);
		integer numslopes;
		if (Nbrstype == NUMBIG)
			numslopes = NNBRS; // h is 6
		else
			numslopes = 11;// h is 4

		q2L = new double[d * m];
		double *q2_coefs = new double[4 * nn];
		a = 0;
		b = 1.0;
		for (integer i = 0; i < d; i++)
		{
			if (isclosed)
			{
				splinestatus = Spline::SplineUniformPeriodic(q2 + i * n, n, 1.0 / nn, q2_coefs);
			}
			else
			{
				splinestatus = Spline::SplineUniformSlopes(q2 + i * n, n, 1.0 / nn, q2_coefs);
			}
			if (!splinestatus)
			{
				printf("Error in computing spline!\n");
				exit(EXIT_FAILURE);
			}
			for (integer j = 0; j < m; j++)
			{
				q2L[j + i * m] = Spline::ValSplineUniform(q2_coefs, n, 1.0 / nn, static_cast<double> (j) / mm);
			}
		}

		delete[] q2_coefs;
		E = new double[n * n];
		for (integer i = 0; i < n * n; i++)
			E[i] = 0;
		Path = new integer[2 * n * n];

		for (integer i = 0; i < n; i++)
		{
			E[n * i + 0] = 1;
			E[n * 0 + i] = 1;
			Path[n * (n * 0 + i) + 0] = -1;
			Path[n * (n * 0 + 0) + i] = -1;
			Path[n * (n * 1 + i) + 0] = -1;
			Path[n * (n * 1 + 0) + i] = -1;
		}
		E[n * 0 + 0] = 0;

		for (integer j = 1; j < n; j++)
		{
			for (integer i = 1; i < n; i++)
			{
				if (Nbrstype == NUMBIG || (j <= 4 * i && i <= 4 * j && i - j < 0.4 * n && j - i < 0.4 * n)) // h is 4
				{
					Emin = 100000;
					Eidx = 0;

					for (Num = 0; Num < numslopes; Num++)
					{
						k = i - Nbrs[Num][0];
						l = j - Nbrs[Num][1];

						if (k >= 0 && l >= 0)
						{
							double slope = static_cast<double> (j - l) / (i - k), sqrts = sqrt(slope), En = 0, y, tmp, tmp2;
							integer idx;
							double dl = static_cast<double> (l);
							integer xmk;

							for (integer x = k; x <= i; x++)
							{
								xmk = x - k;
								y = slope * xmk + dl;
								idx = static_cast<integer> (floor(y * mdn + 0.5));

								tmp2 = 0;
								for (integer h = 0; h < d; h++)
								{
									tmp = q1[x + h * n] - sqrts * q2L[idx + h * m];
									tmp2 += tmp*tmp;
								}
								En += tmp2;
							}

							Etmp = E[n*l + k] + En / (n - 1);

							if (Num == 0 || Etmp < Emin)
							{
								Emin = Etmp;
								Eidx = Num;
							}
						}
					}

					E[n * j + i] = Emin;
					Path[n * (n * 0 + j) + i] = i - Nbrs[Eidx][0];
					Path[n * (n * 1 + j) + i] = j - Nbrs[Eidx][1];
				}
				else
				{
					E[n * j + i] = 100000;
				}
			}
		}
		double Eresult = E[n * n - 1];


		delete[] q2L;
		delete[] E;

		x = new integer[2 * n];
		y = x + n;

		x[0] = n - 1;
		y[0] = n - 1;

		cnt = 1;
		while (x[cnt - 1] > 0)
		{
			y[cnt] = Path[n*(n * 0 + x[cnt - 1]) + y[cnt - 1]];
			x[cnt] = Path[n*(n * 1 + x[cnt - 1]) + y[cnt - 1]];
			cnt++;
		}

		delete[] Path;

		for (integer i = 0, j = cnt - 1; i < j; ++i, --j)
		{
			k = x[i];
			x[i] = x[j];
			x[j] = k;

			k = y[i];
			y[i] = y[j];
			y[j] = k;
		}

		for (integer i = 0; i < n; i++)
		{
			Fmin = 100000;
			Fidx = 0;

			for (integer j = 0; j < cnt; j++)
			{
				Ftmp = (i > x[j] ? i - x[j] : x[j] - i);
				if (j == 0 || Ftmp < Fmin)
				{
					Fmin = Ftmp;
					Fidx = j;
				}
			}

			if (x[Fidx] == i)
			{
				gamma[i] = (y[Fidx]);
			}
			else
			{
				if (x[Fidx] > i)
				{
					a = x[Fidx] - i;
					b = i - x[Fidx - 1];
					gamma[i] = (a * y[Fidx - 1] + b * y[Fidx]) / (a + b);
				}
				else
				{
					a = i - x[Fidx];
					b = x[Fidx + 1] - i;
					gamma[i] = (a * y[Fidx + 1] + b * y[Fidx]) / (a + b);
				}
			}
			gamma[i] /= nn;
		}
		delete[] x;

		return Eresult;
	};

	void CenterC(double *C, integer d, integer n)
	{
		double mean;

		for (integer i = 0; i < d; i++)
		{
			mean = 0;
			for (integer j = 0; j < n; j++)
			{
				mean += C[j + i * n];
			}
			mean /= n;
			for (integer j = 0; j < n; j++)
			{
				C[j + i * n] -= mean;
			}
		}
	};

	void NormalizedC(double *C, integer d, integer n)
	{
		double norm = 0;
		for (integer i = 0; i < d * n; i++)
			norm += C[i] * C[i];
		norm = sqrt(norm);
		for (integer i = 0; i < d * n; i++)
			C[i] /= norm;
	};

	double ComputeTotalAngle(const double *C, integer d, integer n)
	{
		double temp1, temp2, temp3, angle, total_theta = 0;
		double *dx1 = new double[2 * d];
		double *dx2 = dx1 + d;

		for (integer i = 1; i < n - 1; i++)
		{
			temp1 = 0;
			temp2 = 0;
			for (integer j = 0; j < d; j++)
			{
				dx1[j] = C[i + j * n] - C[i - 1 + j * n];
				temp1 += dx1[j] * dx1[j];
				dx2[j] = C[i + 1 + j * n] - C[i + j * n];
				temp2 += dx2[j] * dx2[j];
			}
			temp1 = sqrt(temp1);
			temp2 = sqrt(temp2);
			temp3 = 0;
			for (integer j = 0; j < d; j++)
			{
				dx1[j] /= temp1;
				dx2[j] /= temp2;
				temp3 += dx1[j] * dx2[j];
			}
			angle = acos(temp3);
			if (angle == angle) /*if angle is not nan*/
				total_theta += angle;
		}
		delete[] dx1;
		return total_theta;
	};

	void FindInitialBreaksAndNs(const double *C, integer d, integer n, integer minSkip, double thresholdsmall,
		integer rand_shift, integer *p_ms, integer &Lms, integer &Ns)
	{
		double theta, temp1, temp2, temp3, total_theta, angle, ratio;
		int sind;
		double *dx1 = new double[2 * d + 2 * d * n];
		double *dx2 = dx1 + d;
		double *temppt = dx2 + d;
		for (integer i = 0; i < n - 1; i++)
		{
			for (integer j = 0; j < d; j++)
			{
				temppt[i + j * 2 * n] = C[i + j * n];
				temppt[i + n + j * 2 * n] = C[i + j * n + 1];
			}
		}
		for (integer j = 0; j < d; j++)
		{
			temppt[n - 1 + j * 2 * n] = C[n - 1 + j * n];
		}

		theta = 0;
		total_theta = 0;
		sind = rand_shift;
		(Lms) = 1;
		p_ms[0] = rand_shift;
		ratio = 0;
		for (integer i = 1; i < n - 1; i++)
		{
			//        printf("%f, %f\n", temppt[i + 0 * 2 * n + rand_shift], temppt[i + 1 * 2 * n + rand_shift]);//---
			temp1 = 0;
			temp2 = 0;
			for (integer j = 0; j < d; j++)
			{
				dx1[j] = temppt[i + j * 2 * n + rand_shift] - temppt[i - 1 + j * 2 * n + rand_shift];
				temp1 += dx1[j] * dx1[j];
				dx2[j] = temppt[i + 1 + j * 2 * n + rand_shift] - temppt[i + j * 2 * n + rand_shift];
				temp2 += dx2[j] * dx2[j];
			}
			temp1 = sqrt(temp1);
			temp2 = sqrt(temp2);
			temp3 = 0;
			for (integer j = 0; j < d; j++)
			{
				dx1[j] /= temp1;
				dx2[j] /= temp2;
				temp3 += dx1[j] * dx2[j];
			}
			angle = acos(temp3);
			if (angle == angle) /*if angle is not nan*/
			{
				theta += angle;
				total_theta += angle;
			}
			//        printf("angle: %f, theta: %f\n", angle, theta);//---
			/*        if((dx1 * dy2 - dx2 * dy1 < - sin_angle && dir > 0) || (dx1 * dy2 - dx2 * dy1 > sin_angle && dir < 0))
			threshold = thresholdsmall;*/
			/*        printf("%f, dir: %d, theta: %f, threshold: %f\n", dx1 * dy2 - dx2 * dy1, dir, theta, threshold);*/
			if (theta > thresholdsmall && i + rand_shift - sind >= minSkip)
			{
				theta = 0;
				sind = i + rand_shift;
				/*            p_ms[(*Lms)] = (i + rs >= N) ? i + rs - N : i + rs;*/
				p_ms[Lms] = i + rand_shift;
				Lms++;
			}
		}
		printf("total theta:%g, %f PI\n", total_theta, total_theta / PI);
		Ns = static_cast<int> (static_cast<double> (n) / 3);
		Ns = (Ns > 30) ? 30 : Ns;
		Ns += static_cast<int> (total_theta / PI * 2.0);

		delete[] dx1;
	};

	void CurveToQ(const double *C, integer d, integer n, double *q, bool isclosed)
	{
		double *Ccoefs, *dCcoefs;
		double temp, temp2, tol = sqrt(std::numeric_limits<double>::epsilon());

		Ccoefs = new double[4 * (n - 1) * d + 3 * (n - 1) * d];
		dCcoefs = Ccoefs + 4 * (n - 1) * d;

		for (integer i = 0; i < d; i++)
		{
			if (isclosed)
			{
				Spline::SplineUniformPeriodic(C + i * n, n, 1.0 / (n - 1), Ccoefs + i * 4 * (n - 1));
			}
			else
			{
				Spline::SplineUniformSlopes(C + i * n, n, 1.0 / (n - 1), Ccoefs + i * 4 * (n - 1));
			}
			Spline::FirstDeri(Ccoefs + i * 4 * (n - 1), n, dCcoefs + i * 3 * (n - 1));
		}

		for (integer i = 0; i < n; i++)
		{
			temp = 0;
			for (integer j = 0; j < d; j++)
			{
				q[i + j * n] = Spline::ValFirstDeriUniform(dCcoefs + j * 3 * (n - 1), n, 1.0 / (n - 1), static_cast<double> (i) / (n - 1));
				temp += q[i + j * n] * q[i + j * n];
			}
			temp = sqrt(temp);

			if (temp > tol)
			{
				for (integer j = 0; j < d; j++)
				{
					q[i + j * n] /= temp;
				}
			}
			else
			{
				for (integer j = 0; j < d; j++)
				{
					q[i + j * n] = 0;
				}
			}
		}
		temp = 0;
		for (integer i = 0; i < n; i++)
		{
			temp2 = 0;
			for (integer j = 0; j < d; j++)
			{
				temp2 += q[i + j * n] * q[i + j * n];
			}
			if (i == 0 || i == n - 1)
			{
				temp2 /= 2;
			}
			temp += temp2;
		}
		temp /= n - 1;
		temp = sqrt(temp);
		for (integer i = 0; i < d * n; i++)
		{
			q[i] /= temp;
		}
		delete[] Ccoefs;
	};


	void QToCurve(const double *Q, integer d, integer n, double *C, bool isclosed)
	{
		double *q2n = new double[n + n * d];
		double *q2nTimesQ = q2n + n;
		integer inc = n;

		for (integer i = 0; i < n; i++)
		{
			q2n[i] = sqrt(ddot_(&d, const_cast<double *> (Q + i), &inc, const_cast<double *> (Q + i), &inc));
		}

		ElasticCurvesRO::PointwiseQProdl(Q, q2n, d, n, q2nTimesQ);
		for (integer i = 0; i < d; i++)
		{
			ElasticCurvesRO::CumTrapz(q2nTimesQ + i * n, n, 1.0 / (n - 1), C + i * n);
		}

		delete[] q2n;
	};

	void ShiftC(const double *C, integer d, integer n, double *Cshift, integer m)
	{
		if (m == 0) // in case of open curve
		{
			integer nd = n * d, inc = 1;
			dcopy_(&nd, const_cast<double *> (C), &inc, Cshift, &inc);
			return;
		}

		int ind, Nm1;
		Nm1 = n - 1;
		for (integer i = 0; i < Nm1; i++)
		{
			ind = i + Nm1 - m;
			ind = (ind >= Nm1) ? ind - Nm1 : ind;
			for (integer j = 0; j < d; j++)
			{
				Cshift[ind + j * n] = C[i + j * n];
			}
		}
		for (integer j = 0; j < d; j++)
		{
			Cshift[Nm1 + j * n] = Cshift[j * n];
		}
	};

	void FindBestRotation(const double *q1, const double *q2, integer d, integer n, double *O)
	{
		double *M = new double[d * d]; // M = q1^T G q2
		for (integer i = 0; i < d; i++)
		{
			for (integer j = 0; j < d; j++)
			{
				M[i + j * d] = q1[i * n] * q2[j * n] / 2;
				for (integer k = 1; k < n - 1; k++)
				{
					M[i + j * d] += q1[k + i * n] * q2[k + j * n];
				}
				M[i + j * d] += q1[n - 1 + i * n] * q2[n - 1 + j * n] / 2;
				M[i + j * d] /= (n - 1);
			}
		}

		// compute SVD of M;
		char *joba = const_cast<char *> ("A");
		double *U = new double[2 * d * d + d];
		double *Vt = U + d * d;
		double *S = Vt + d * d;
		integer lwork = -1, info;
		double workoptsize;
		// compute the space required in the SVD computation.
		dgesvd_(joba, joba, &d, &d, M, &d, S, U, &d, Vt, &d, &workoptsize, &lwork, &info);
		lwork = static_cast<integer> (workoptsize);
		double *work = new double[lwork];
		// SVD: U * S * Vt = M, details: http://www.netlib.org/lapack/explore-html/d8/d2d/dgesvd_8f.html
		dgesvd_(joba, joba, &d, &d, M, &d, S, U, &d, Vt, &d, work, &lwork, &info);
		if (info != 0)
		{
			printf("Error:singular value decomposition failed!\n");
		}
		delete[] M;
		delete[] work;
		char *transn = const_cast<char *> ("n");
		double one = 1, zero = 0;
		// O = U * Vt, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &d, &d, &d, &one, U, &d, Vt, &d, &zero, O, &d);

		// obtain the determinant of O
		double *O2 = new double[d * d];
		for (integer i = 0; i < d * d; i++)
			O2[i] = O[i];

		//===========evaluate det(O2)================
		integer *IPIV = new integer[d];
		// LU decomposition for O2, P * L * U <- A, L and U are stored in A, permutation matrix P are in IPIV
		// details:http://www.netlib.org/lapack/explore-html/d3/d6a/dgetrf_8f.html
		dgetrf_(&d, &d, O2, &d, IPIV, &info);
		if (info != 0)
		{
			printf("Error:LU decomposition failed!\n");
		}

		double det = 1;
		for (integer i = 0; i < d; i++)
		{
			if (IPIV[i] != i + 1)
			{
				det *= -O2[i + i * d];
			}
			else
			{
				det *= O2[i + i * d];
			}
		}
		delete[] IPIV;
		delete[] O2;
		//===========finish evaluating det(O2)================


		if (det > 0)
		{
			delete[] U;
			return;
		}

		for (integer i = 0; i < d; i++)
		{
			Vt[d - 1 + i * d] *= -1.0;
		}

		// O = U * Vt, details: http://www.netlib.org/lapack/explore-html/d7/d2b/dgemm_8f.html
		dgemm_(transn, transn, &d, &d, &d, &one, U, &d, Vt, &d, &zero, O, &d);

		delete[] U;
	};

	void GetCurveSmall(const double *C, double *Cs, integer d, integer n, integer ns, bool isclosed)
	{
		double *coefs;
		coefs = new double[4 * (n - 1)];

		for (integer i = 0; i < d; i++)
		{
			if (isclosed)
			{
				Spline::SplineUniformPeriodic(C + i * n, n, 1.0 / (n - 1), coefs);
			}
			else
			{
				Spline::SplineUniformSlopes(C + i * n, n, 1.0 / (n - 1), coefs);
			}
			for (integer j = 0; j < ns; j++)
			{
				Cs[j + i * ns] = Spline::ValSplineUniform(coefs, n, 1.0 / (n - 1), static_cast<double> (j) / (ns - 1));
			}
		}

		delete[] coefs;
	};

	void GammaInverse(const double *DPgam, integer n, double *DPgamI)
	{
		double *DPIcoefs, *breaks;
		int Nm1, ind;
		Nm1 = n - 1;
		DPIcoefs = new double[4 * Nm1 + n];
		breaks = DPIcoefs + 4 * Nm1;
		for (integer i = 0; i < n; i++)
			breaks[i] = static_cast<double> (i) / Nm1;

		DPgamI[0] = 0;
		ind = 0;
		for (integer i = 1; i < n; i++)
		{
			while (breaks[i] > DPgam[ind + 1] && ind + 1 < n)
				ind++;
			DPgamI[i] = breaks[ind] + (breaks[ind + 1] - breaks[ind]) / (DPgam[ind + 1] - DPgam[ind]) * (breaks[i] - DPgam[ind]);
		}
		delete[] DPIcoefs;
	};

	void ReSampleGamma(const double *DPgams, integer ns, double *DPgam, integer n)
	{
		double *coefs, *breaks;
		double x;
		int ind;
		coefs = new double[4 * (ns - 1) + ns];
		breaks = coefs + 4 * (ns - 1);

		for (integer i = 0; i < ns; i++)
			breaks[i] = static_cast<double> (i) / (ns - 1);

		DPgam[0] = DPgams[0];
		ind = 0;
		for (integer i = 1; i < n; i++)
		{
			x = static_cast<double> (i) / (n - 1);
			while (x > breaks[ind + 1] && ind + 1 < ns)
				ind++;
			DPgam[i] = DPgams[ind] + (DPgams[ind + 1] - DPgams[ind]) / (breaks[ind + 1] - breaks[ind]) * (x - breaks[ind]);
		}

		delete[] coefs;
	};

	void Gradient(const double *DPgam, integer n, double h, double *grad)
	{
		grad[0] = (DPgam[1] - DPgam[0]) / h;
		grad[n - 1] = (DPgam[n - 1] - DPgam[n - 2]) / h;

		for (integer i = 1; i < n - 1; i++)
			grad[i] = (DPgam[i + 1] - DPgam[i - 1]) / 2 / h;
	};

	void GradientPeriod(const double *DPgam, integer n, double h, double *grad)
	{
		int i;
		grad[0] = (DPgam[1] - DPgam[0] + DPgam[n - 1] - DPgam[n - 2]) / 2 / h;
		grad[n - 1] = grad[0];
		for (i = 1; i < n - 1; i++)
			grad[i] = (DPgam[i + 1] - DPgam[i - 1]) / 2 / h;
	};
}; /*end of ROPTLIB namespace*/
