
#include "Others/Spline.h"

/*Define the namespace*/
namespace ROPTLIB{

	int Spline::SplineUniformPeriodic(const double *Y, int n, double h, double *coefs)
	{ // solving system based on second derivatives.
		int i, nn;
		double *d, *ud, *ld, *vec, *s;
		nn = n - 1;
		d = new double[5 * nn - 1];
		ud = d + nn;
		ld = ud + (nn - 1);
		vec = ld + (nn - 1);
		s = vec + nn;
		if (fabs(Y[0] - Y[nn]) > sqrt(std::numeric_limits<double>::epsilon()))
		{
			printf("warning: %g = Y[start] != Y[end] = %g: %g, Using curbic spline with periodic condition may cause problems.\n", Y[0], Y[nn], Y[0] - Y[nn]);
		}

		for (i = 0; i < nn; i++)
		{
			ld[i] = 0.5;
			d[i] = 2.0;
			ud[i] = 0.5;
			if (i == nn - 1)
				vec[i] = 3.0 / h * ((Y[1] - Y[i + 1]) / h - (Y[i + 1] - Y[i]) / h);
			else
				vec[i] = 3.0 / h * ((Y[i + 2] - Y[i + 1]) / h - (Y[i + 1] - Y[i]) / h);
		}

		if (!SolvePeriodicSystem(d, ud, ld, vec, s, nn))
		{
			printf("error: fail to slove the linear system!!\n");
			return 0;
		}

		s[0] = s[nn];
		for (i = 0; i < nn; i++)
		{
			coefs[0 * nn + i] = (s[i + 1] - s[i]) / 6 / h;
			coefs[1 * nn + i] = s[i] / 2;
			coefs[2 * nn + i] = -s[i] * h / 2 + (Y[i + 1] - Y[i]) / h - h * (s[i + 1] - s[i]) / 6;
			coefs[3 * nn + i] = Y[i];
		}
		delete[] d;
		return 1;
	};

	int Spline::SplinePeriodic(const double *X, const double *Y, int n, double *coefs)
	{ // solving system based on second derivatives.
		int i, nn;
		double *d, *ud, *ld, *vec, *s;
		double hi;
		nn = n - 1;
		d = new double[5 * nn - 1];
		ud = d + nn;
		ld = ud + (nn - 1);
		vec = ld + (nn - 1);
		s = vec + nn;
		if (fabs(Y[0] - Y[nn]) > sqrt(std::numeric_limits<double>::epsilon()))
		{
			printf("warning: %g = Y[start] != Y[end] = %g, %g, Using curbic spline with periodic condition may cause problems.\n", Y[0], Y[nn], Y[0] - Y[nn]);
		}
		for (i = 0; i < nn; i++)
		{
			if (i == nn - 1)
			{
				ld[i] = (X[i + 1] - X[i]) / (X[1] - X[0] + X[i + 1] - X[i]);
				d[i] = 2;
				ud[i] = (X[1] - X[0]) / (X[1] - X[0] + X[i + 1] - X[i]);
				vec[i] = 6 / ((X[1] - X[0] + X[i + 1] - X[i])) * ((Y[1] - Y[0]) / (X[1] - X[0]) - (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]));
			}
			else
			{
				ld[i] = (X[i + 1] - X[i]) / (X[i + 2] - X[i]);
				d[i] = 2;
				ud[i] = (X[i + 2] - X[i + 1]) / (X[i + 2] - X[i]);
				vec[i] = 6 / (X[i + 2] - X[i]) * ((Y[i + 2] - Y[i + 1]) / (X[i + 2] - X[i + 1]) - (Y[i + 1] - Y[i]) / (X[i + 1] - X[i]));
			}
		}

		if (!SolvePeriodicSystem(d, ud, ld, vec, s, nn))
		{
			printf("error: fail to slove the linear system!!\n");
			return 0;
		}

		s[0] = s[nn];
		for (i = 0; i < nn; i++)
		{
			hi = (X[i + 1] - X[i]);
			coefs[0 * nn + i] = (s[i + 1] - s[i]) / 6 / hi;
			coefs[1 * nn + i] = s[i] / 2;
			coefs[2 * nn + i] = -s[i] * hi / 2 + (Y[i + 1] - Y[i]) / hi - hi * (s[i + 1] - s[i]) / 6;
			coefs[3 * nn + i] = s[i] * hi * hi / 6 + (Y[i] - s[i] * hi * hi / 6);
		}

		delete[] d;
		return 1;
	};

	int Spline::SplineUniformSlopes(const double *Y, int n, double h, double *coefs)
	{ // solving system based on second derivatives.
		int i, nn;
		double *d, *ud, *ld, *vec, *s;
		double first_d;
		d = new double[5 * n - 2];
		ud = d + n;
		ld = ud + (n - 1);
		vec = ld + (n - 1);
		s = vec + n;

		nn = n - 1;
		for (i = 1; i < nn; i++)
		{
			ld[i - 1] = 0.5;
			d[i] = 2.0;
			ud[i] = 0.5;
			vec[i] = 3.0 / h * ((Y[i + 1] - Y[i]) / h - (Y[i] - Y[i - 1]) / h);
		}
		/*s2 form*/
		first_d = (Y[1] - Y[0]) / (h);
		d[0] = h / 3;
		ud[0] = h / 6;
		vec[0] = (Y[1] - Y[0]) / h - first_d;

		first_d = (Y[nn] - Y[nn - 1]) / (h);
		d[nn] = h / 3;
		ld[nn - 1] = h / 6;
		vec[nn] = first_d - (Y[nn] - Y[nn - 1]) / h;

		if (!SolveTridiagonalSystem(d, ud, ld, vec, s, n))
		{
			printf("error: fail to slove tridiagonal system!!\n");
			return 0;
		}

		for (i = 0; i < nn; i++)
		{
			coefs[0 * nn + i] = (s[i + 1] - s[i]) / 6 / h;
			coefs[1 * nn + i] = s[i] / 2;
			coefs[2 * nn + i] = -s[i] * h / 2 + (Y[i + 1] - Y[i]) / h - h * (s[i + 1] - s[i]) / 6;
			coefs[3 * nn + i] = Y[i];
		}
		delete[] d;
		return 1;
	};

	int Spline::SplineSlopes(const double *X, const double *Y, int n, double *coefs)
	{ // solving system based on second derivatives.
		int i, nn;
		double *d, *ud, *ld, *vec, *s;
		double first_d, hi;
		d = new double[5 * n - 2];
		ud = d + n;
		ld = ud + (n - 1);
		vec = ld + (n - 1);
		s = vec + n;
		nn = n - 1;
		for (i = 1; i < nn; i++)
		{
			ld[i - 1] = (X[i] - X[i - 1]) / (X[i + 1] - X[i - 1]);
			d[i] = 2;
			ud[i] = (X[i + 1] - X[i]) / (X[i + 1] - X[i - 1]);
			vec[i] = 6 / (X[i + 1] - X[i - 1]) * ((Y[i + 1] - Y[i]) / (X[i + 1] - X[i]) - (Y[i] - Y[i - 1]) / (X[i] - X[i - 1]));
		}
		/*s2 form*/
		first_d = (Y[1] - Y[0]) / (X[1] - X[0]);
		d[0] = (X[1] - X[0]) / 3;
		ud[0] = (X[1] - X[0]) / 6;
		vec[0] = (Y[1] - Y[0]) / (X[1] - X[0]) - first_d;

		first_d = (Y[nn] - Y[nn - 1]) / (X[nn] - X[nn - 1]);
		d[nn] = (X[nn] - X[nn - 1]) / 3;
		ld[nn - 1] = (X[nn] - X[nn - 1]) / 6;
		vec[nn] = first_d - (Y[nn] - Y[nn - 1]) / (X[nn] - X[nn - 1]);

		if (!SolveTridiagonalSystem(d, ud, ld, vec, s, n))
		{
			printf("error: fail to slove tridiagonal system!!\n");
			return 0;
		}

		for (i = 0; i < nn; i++)
		{
			hi = (X[i + 1] - X[i]);
			coefs[0 * nn + i] = (s[i + 1] - s[i]) / 6 / hi;
			coefs[1 * nn + i] = s[i] / 2;
			coefs[2 * nn + i] = -s[i] * hi / 2 + (Y[i + 1] - Y[i]) / hi - hi * (s[i + 1] - s[i]) / 6;
			coefs[3 * nn + i] = s[i] * hi * hi / 6 + (Y[i] - s[i] * hi * hi / 6);
		}

		delete[] d;
		return 1;
	};

	int Spline::SolveTridiagonalSystem(double *d, double *ud, double *ld, double *vec, double *s, int n)
	{
		int i;
		double coef;
		for (i = 0; i < n - 1; i++)
		{
			coef = -ld[i] / d[i];
			ld[i] += d[i] * coef;
			d[i + 1] += ud[i] * coef;
			vec[i + 1] += vec[i] * coef;
		}
		if (fabs(d[n - 1]) < std::numeric_limits<double>::epsilon())
		{
			printf("tridiagonal system can not be solved!!");
			return 0;
		}
		s[n - 1] = vec[n - 1] / d[n - 1];
		for (i = n - 2; i >= 0; i--)
		{
			if (fabs(d[i]) < std::numeric_limits<double>::epsilon())
			{
				printf("tridiagonal system can not be solved!!");
				return 0;
			}
			s[i] = (vec[i] - s[i + 1] * ud[i]) / d[i];
		}
		return 1;
	};

	int Spline::SolvePeriodicSystem(double *d, double *ud, double *ld, double *vec, double *s, int nn)
	{
		double temp = ud[nn - 1];
		double coef;
		int i;
		double *last_column = new double[nn - 2];
		last_column[0] = ld[0];

		for (i = 0; i < nn - 3; i++)
		{
			coef = -ld[i + 1] / d[i];
			d[i + 1] += ud[i] * coef;
			last_column[i + 1] = last_column[i] * coef;
			vec[i + 1] += vec[i] * coef;

			coef = -temp / d[i]; /*this temp stores nn - 1 row, i column entry*/
			temp = ud[i] * coef; /*this temp stores nn - 1 row, i + 1 column entry*/
			d[nn - 1] += last_column[i] * coef;
			vec[nn - 1] += vec[i] * coef;
		}

		/*i = nn - 3*/
		i = nn - 3;
		coef = -ld[i + 1] / d[i];
		d[i + 1] += ud[i] * coef;
		ud[nn - 2] += last_column[i] * coef;
		vec[i + 1] += vec[i] * coef;

		coef = -temp / d[i]; /*this temp stores nn - 1 row, i column entry*/
		ld[nn - 1] += ud[i] * coef; /*this temp stores nn - 1 row, i + 1 column entry*/
		d[nn - 1] += last_column[i] * coef;
		vec[nn - 1] += vec[i] * coef;

		/*i = nn - 2*/
		i = nn - 2;
		coef = -ld[i + 1] / d[i];
		d[i + 1] += ud[i] * coef;
		vec[i + 1] += vec[i] * coef;

		/*solve upper triangle problem*/
		s[nn] = vec[nn - 1] / d[nn - 1];
		if (fabs(d[nn - 1]) < std::numeric_limits<double>::epsilon())
		{
			printf("upper triangle system can not be solved!!");
			return 0;
		}
		s[nn - 1] = (vec[nn - 2] - s[nn] * ud[nn - 2]) / d[nn - 2];
		for (i = nn - 2; i > 0; i--)
		{
			if (fabs(d[i - 1]) < std::numeric_limits<double>::epsilon())
			{
				printf("upper triangle system can not be solved!!");
				return 0;
			}
			s[i] = (vec[i - 1] - s[nn] * last_column[i - 1] - s[i + 1] * ud[i - 1]) / d[i - 1];
		}
		s[0] = s[nn];

		delete[] last_column;
		return 1;
	};

	double Spline::ValSpline(const double *coefs, const double *breaks, int N, double t)
	{
		int i, nn;
		double output;
		nn = N - 1;
		i = 0;
		while (i < N && t - (breaks[i] - breaks[0]) >= -std::numeric_limits<double>::epsilon())
			i++;
		i--;
		i = (i < 0) ? 0 : i;
		i = (i > nn - 1) ? nn - 1 : i;
		t -= breaks[i];
		output = ((coefs[0 * nn + i] * t + coefs[1 * nn + i]) * t + coefs[2 * nn + i]) * t + coefs[3 * nn + i];
		return output;
	};

	double Spline::ValSplineUniform(const double *coefs, int N, double h, double t)
	{
		integer i, nn;
		double output;
		nn = N - 1;
		i = static_cast<integer> (t / h);
		while (t - i * h >= -std::numeric_limits<double>::epsilon())
			i++;
		i--;
		i = (i < 0) ? 0 : i;
		i = (i > nn - 1) ? nn - 1 : i;
		t -= i * h;
		output = ((coefs[0 * nn + i] * t + coefs[1 * nn + i]) * t + coefs[2 * nn + i]) * t + coefs[3 * nn + i];
		return output;
	};

	void Spline::FirstDeri(const double *coefs, int N, double *dericoefs)
	{
		int i, nn;
		nn = N - 1;
		for (i = 0; i < nn; i++)
		{
			dericoefs[0 * nn + i] = coefs[0 * nn + i] * 3;
			dericoefs[1 * nn + i] = coefs[1 * nn + i] * 2;
			dericoefs[2 * nn + i] = coefs[2 * nn + i];
		}
	};

	double Spline::ValFirstDeriUniform(const double *dericoefs, int N, double h, double t)
	{
		int i, nn;
		nn = N - 1;
		i = static_cast<int> (t / h);
		while (t - i * h >= -std::numeric_limits<double>::epsilon())
			i++;
		i--;
		i = (i < 0) ? 0 : i;
		i = (i > nn - 1) ? nn - 1 : i;
		t -= i * h;
		return (dericoefs[0 * nn + i] * t + dericoefs[1 * nn + i]) * t + dericoefs[2 * nn + i];
	};

	double Spline::ValFirstDeri(const double *dericoefs, const double *breaks, int N, double t)
	{
		int i, nn;
		double output;
		nn = N - 1;
		i = 0;

		while (i < N && t - (breaks[i] - breaks[0]) >= -std::numeric_limits<double>::epsilon())
			i++;
		i--;
		i = (i < 0) ? 0 : i;
		i = (i > nn - 1) ? nn - 1 : i;
		t -= breaks[i];
		output = (dericoefs[0 * nn + i] * t + dericoefs[1 * nn + i]) * t + dericoefs[2 * nn + i];
		return output;
	};

	void Spline::SecondDeri(const double *coefs, int N, double *dericoefs)
	{
		int i, nn;
		nn = N - 1;
		for (i = 0; i < nn; i++)
		{
			dericoefs[0 * nn + i] = coefs[0 * nn + i] * 6;
			dericoefs[1 * nn + i] = coefs[1 * nn + i] * 2;
		}
	};

	double Spline::ValSecondDeriUniform(const double *dericoefs, int N, double h, double t)
	{
		int i, nn;
		nn = N - 1;
		i = static_cast<int> (t / h);
		while (t - i * h >= -std::numeric_limits<double>::epsilon())
			i++;
		i--;
		i = (i < 0) ? 0 : i;
		i = (i > nn - 1) ? nn - 1 : i;
		t -= i * h;
		return dericoefs[0 * nn + i] * t + dericoefs[1 * nn + i];
	};

	double Spline::ValSecondDeri(const double *dericoefs, const double *breaks, int N, double t)
	{
		int i, nn;
		double output;
		nn = N - 1;
		i = 0;
		while (i < N && t - (breaks[i] - breaks[0]) >= -std::numeric_limits<double>::epsilon())
			i++;
		i--;
		i = (i < 0) ? 0 : i;
		i = (i > nn - 1) ? nn - 1 : i;
		t -= breaks[i];
		output = dericoefs[0 * nn + i] * t + dericoefs[1 * nn + i];
		return output;
	};
}; /*end of ROPTLIB namespace*/
