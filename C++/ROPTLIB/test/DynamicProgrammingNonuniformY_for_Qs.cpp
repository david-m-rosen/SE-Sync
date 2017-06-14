
#include "DynamicProgrammingNonuniformY.h"

#ifdef MATLAB_MEX_FILE

std::map<integer *, integer> *CheckMemoryDeleted;

double DynamicProgrammingNonuniformY(const double *q1, const double *q2, const double *q2y, integer d, integer n, double *gamma, bool isclosed);
void CurveToQNonUniformY(const double *Cgrid, const double *C, integer d, integer n, double *q, bool isclosed);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 4)
	{
		mexErrMsgTxt("The number of arguments should be 4.\n");
	}
	double *q1, *q2, *q2grid;
    integer isclosed;
//     printf("s0\n");//---
	q1 = mxGetPr(prhs[0]);
	q2 = mxGetPr(prhs[1]);
    q2grid = mxGetPr(prhs[2]);
	isclosed = static_cast<integer> (mxGetScalar(prhs[3]));
	/* dimensions of input matrices */
	integer d, n;
	n = mxGetM(prhs[0]);
	d = mxGetN(prhs[0]);
//     printf("s1\n");//---

	std::cout << "(n, d):" << n << "," << d << std::endl;

	if (mxGetM(prhs[1]) != n || mxGetN(prhs[1]) != d)
	{
		mexErrMsgTxt("The size of matrix q2 does not match the size of q1.\n");
	}
    if(mxGetM(prhs[2]) != n)
    {
		mexErrMsgTxt("The size of matrix q2grid does not match the size of q2.\n");
    }
    
//     printf("s2\n");//---
	genrandseed(0);

//     printf("s21\n");//---
	CheckMemoryDeleted = new std::map<integer *, integer>;

//     printf("s22\n");//---
    
//     double *q1, *q2;
// //     printf("s23\n");//---
//     q1 = new double [2 * d * n];
// //     printf("s24\n");//---
//     q2 = q1 + d * n;
// //     printf("s25\n");//---
// 	CurveToQ(C1, d, n, q1, isclosed == 1);
// //     printf("s26\n");//---
// 	CurveToQNonUniformY(C2grid, C2, d, n, q2, isclosed == 1);
//     
//     printf("s3\n");//---
	/*create output matrix*/
	plhs[0] = mxCreateDoubleMatrix(n, 1, mxREAL);
	double *optgam = mxGetPr(plhs[0]);
    
//     printf("s4\n");//---
    double energy = DynamicProgrammingNonuniformY(q1, q2, q2grid, d, n, optgam, isclosed == 1);
	plhs[1] = mxCreateDoubleScalar(energy);
    
	std::map<integer *, integer>::iterator iter = CheckMemoryDeleted->begin();
	for (iter = CheckMemoryDeleted->begin(); iter != CheckMemoryDeleted->end(); iter++)
	{
		if (iter->second != 1)
			std::cout << "Global address:" << iter->first << ", sharedtimes:" << iter->second << std::endl;
	}
	delete CheckMemoryDeleted;
//     delete [] q1;
	return;
};

double DynamicProgrammingNonuniformY(const double *q1, const double *q2, const double *q2y, integer d, integer n, double *gamma, bool isclosed)
{
	integer k, l, m, Eidx, Fidx, Ftmp, Fmin, Num, *Path, *x, *y, cnt;
	double *q2L, *E, Etmp, Emin, a, b;
	integer nn = n - 1, splinestatus;
	m = 5 * nn + 1;
	integer mm = m - 1;
//     printf("t0\n");//---
	double mdn = static_cast<double> (mm) / (nn);
	q2L = new double[d * m];
	double *q2_coefs = new double[4 * nn];
	a = 0;
	b = 1.0;
//     printf("t1\n");//---
	for (integer i = 0; i < d; i++)
	{
		if (isclosed)
		{
            splinestatus = Spline::SplineUniformPeriodic(q2 + i * n, n, 1.0 / nn, q2_coefs);//--(const double *Y, int n, double h, double *coefs)
		}
		else
		{
            splinestatus = Spline::SplineUniformSlopes(q2 + i * n, n, 1.0 / nn, q2_coefs);
		}
		if (!splinestatus)
		{
			std::cout << "Error in computing spline!" << std::endl;
			exit(EXIT_FAILURE);
		}
        
        integer idx = 1;
        double gamION = 0;
		for (integer j = 0; j < m; j++)
		{
            while(idx < n && static_cast<double> (j) / mm > q2y[idx])
            {
                idx++;
            }
            gamION = static_cast<double> (idx - 1) / nn + 1.0 / nn / (q2y[idx] - q2y[idx - 1]) * (static_cast<double> (j) / mm - q2y[idx - 1]);
//             printf("idx: %d, j: %d, gamION: %f, nn:%d, mm:%d , q2y[idx]:%f, q2y[idx-1]:%f\n", idx, j, gamION, nn, mm, q2y[idx], q2y[idx-1]);//----
            q2L[j + i * m] =  Spline::ValSplineUniform(q2_coefs, n, 1.0 / nn, gamION);
		}
	}

//     printf("t2\n");//---
	delete[] q2_coefs;
	E = new double[n * n];
	for (integer i = 0; i < n * n; i++)
		E[i] = 0;
	Path = new integer[2 * n * n];

//     printf("t3\n");//---
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

//     printf("t4\n");//---
	for (integer j = 1; j < n; j++)
	{
		for (integer i = 1; i < n; i++)
		{
			Emin = 100000;
			Eidx = 0;

			for (Num = 0; Num < NNBRSNON; Num++)
			{
				k = i - NbrsNon[Num][0];
				l = j - NbrsNon[Num][1];

				if (k >= 0 && l >= 0)
				{
					double slope = static_cast<double> (q2y[j] - q2y[l]) * n / (i - k), sqrts = sqrt(slope), En = 0, y, tmp, tmp2;
					integer idx;
// 					double dl = static_cast<double> (l);
					double xmk;

					for (integer x = k; x <= i; x++)
					{
						xmk = static_cast<double> (x - k) / n;
						y = slope * xmk + q2y[l];
						idx = static_cast<integer> (floor(y * mm + 0.5));

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
			Path[n * (n * 0 + j) + i] = i - NbrsNon[Eidx][0];
			Path[n * (n * 1 + j) + i] = j - NbrsNon[Eidx][1];
		}
	}
//     printf("t5\n");//---
	double Eresult = E[n * n - 1];


//     printf("t6\n");//---
	delete[] q2L;
//     printf("t7\n");//---
	delete[] E;

	x = new integer[2 * n];
	y = x + n;

	x[0] = n - 1;
	y[0] = n - 1;

//     printf("t8\n");//---
	cnt = 1;
	while (x[cnt - 1] > 0)
	{
		x[cnt] = Path[n*(n * 0 + y[cnt - 1]) + x[cnt - 1]];
		y[cnt] = Path[n*(n * 1 + y[cnt - 1]) + x[cnt - 1]];
		cnt++;
	}

//     printf("t9\n");//---
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

//     printf("t10\n");//---
//     printf("cnt:%d\n", cnt);//---
//     printf("n:%d\n", n);//---
    double *yd = new double [cnt];
    
    for (integer i = 0; i < cnt; i++)
    {
//         printf("i:%d, %d\n", i, y[i]);//---
        yd[i] = q2y[static_cast<integer> (y[i])];
//         printf("i:%d, %f\n", i, yd[i]);//---
    }
    
//     printf("t11\n");//---
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
			gamma[i] = (yd[Fidx]);
		}
		else
		{
			if (x[Fidx] > i)
			{
				a = x[Fidx] - i;
				b = i - x[Fidx - 1];
				gamma[i] = (a * yd[Fidx - 1] + b * yd[Fidx]) / (a + b);
			}
			else
			{
				a = i - x[Fidx];
				b = x[Fidx + 1] - i;
				gamma[i] = (a * yd[Fidx + 1] + b * yd[Fidx]) / (a + b);
			}
		}
//         gamma[i] /= nn;//--
	}
//     printf("t12\n");//---
	delete[] x;
    delete [] yd;
    
//     printf("t13\n");//---
//     for (integer i = 0; i < n; i++)
//     {
//         printf("i:%d, %f\n", i, gamma[i] * nn);//----
//     }

	return Eresult;
};

void CurveToQNonUniformY(const double *Cgrid, const double *C, integer d, integer n, double *q, bool isclosed)
{
	double *Ccoefs, *dCcoefs;
	double temp, temp1, temp2, tol = sqrt(std::numeric_limits<double>::epsilon());

//     printf("q1\n");//---
	Ccoefs = new double[4 * (n - 1) * d + 3 * (n - 1) * d];
	dCcoefs = Ccoefs + 4 * (n - 1) * d;

//     printf("q2\n");//---
	for (integer i = 0; i < d; i++)
	{
		if (isclosed)
		{
			Spline::SplinePeriodic(Cgrid, C + i * n, n, Ccoefs + i * 4 * (n - 1));
		}
		else
		{
			Spline::SplineSlopes(Cgrid, C + i * n, n, Ccoefs + i * 4 * (n - 1));
		}
		Spline::FirstDeri(Ccoefs + i * 4 * (n - 1), n, dCcoefs + i * 3 * (n - 1));
	}
//     printf("q3\n");//---

	for (integer i = 0; i < n; i++)
	{
		temp = 0;
//     printf("q31\n");//---
		for (integer j = 0; j < d; j++)
		{
			q[i + j * n] = Spline::ValFirstDeri(dCcoefs + j * 3 * (n - 1), Cgrid, n, Cgrid[i]);
			temp += q[i + j * n] * q[i + j * n];
		}
		temp = sqrt(temp);

//     printf("q32\n");//---
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
//     printf("q4\n");//---
	temp = 0;
	for (integer i = 0; i < n - 1; i++)
	{
        temp1 = 0;
		temp2 = 0;
		for (integer j = 0; j < d; j++)
		{
			temp1 += q[i + j * n] * q[i + j * n];
            temp2 += q[i + 1 + j * n] * q[i + 1 + j * n];
		}
        temp += (temp1 + temp2) * (Cgrid[i + 1] - Cgrid[i]) / 2;
	}
//     printf("q5\n");//---
	temp = sqrt(temp);
	for (integer i = 0; i < d * n; i++)
	{
		q[i] /= temp;
	}
//     printf("q6\n");//---
	delete[] Ccoefs;
//     printf("q7\n");//---
};

#endif
