/*
Some blas and lapack functions are not easy to use. This class defines wrapper functions
which can be used to solve some problems or decomposition, e.g., SVD decomposition, Sylevster equation.
More functions will be added in this class.

-----WH
*/

#ifndef MYMATRIX_H
#define MYMATRIX_H

//#include <cmath>
#include <ctime>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cassert>
#include "Others/def.h"

/*Define the namespace*/
namespace ROPTLIB{

	namespace GLOBAL{
		extern integer IZERO, IONE, ITWO;
		extern double DZERO, DONE, DTWO, DNONE;
		extern doublecomplex ZZERO, ZONE, ZTWO, ZNONE;
		extern char *N, *T, *L, *R, *V, *C, *U, *A, *S;
	};

	class Matrix{


		friend std::istream &operator>>(std::istream &input, Matrix &mat)
		{
			return input;
		};

		friend std::ostream &operator<<(std::ostream &output, const Matrix &mat)
		{
			printf("{(%d, %d)\n", mat.row, mat.col);
			for (int i = 0; i < mat.row; i++)
			{
				for (int j = 0; j < mat.col; j++)
					if (j != mat.col - 1)
						printf("%f, ", mat.matrix[i + j * mat.inc]);
					else
						printf("%f", mat.matrix[i + j * mat.inc]);
				if (i != mat.row - 1)
					printf("\n");
			}
			printf("\n}\n");

			return output;
		};

		// All below overload functions are intuitive and easy to use. However, they all have low efficiency.
		// Only use for DEBUG!!=================================================================================

		friend Matrix operator+(const Matrix &left, const Matrix &right) // low efficient
		{
			printf("For debug, low efficient!\n");
			assert(left.row == right.row && right.col == left.col);
			Matrix result(left.row, left.col);

			for (int i = 0; i < left.row; i++)
				for (int j = 0; j < left.col; j++)
					result.matrix[i + j * result.inc] = left.matrix[i + j * left.inc] + right.matrix[i + j * right.inc];

			return result;
		};

		friend Matrix operator+(const Matrix &left, const double &right) // low efficient
		{
			printf("For debug, low efficient!\n");
			Matrix result(left.row, left.col);

			for (int i = 0; i < left.row; i++)
				for (int j = 0; j < left.col; j++)
					result.matrix[i + j * result.inc] = left.matrix[i + j * left.inc] + right;

			return result;
		};

		friend Matrix operator+(const double &left, const Matrix &right)// low efficient
		{
			printf("For debug, low efficient!\n");
			Matrix result(right.row, right.col);

			for (int i = 0; i < right.row; i++)
				for (int j = 0; j < right.col; j++)
					result.matrix[i + j * result.inc] = left + right.matrix[i + j * right.inc];

			return result;
		};

		friend Matrix operator-(const Matrix &left, const Matrix &right)// low efficient
		{
			printf("For debug, low efficient!\n");
			assert(left.row == right.row && right.col == left.col);

			Matrix result(right.row, right.col);

			for (int i = 0; i < left.row; i++)
				for (int j = 0; j < left.col; j++)
					result.matrix[i + j * result.inc] = left.matrix[i + j * left.inc] - right.matrix[i + j * right.inc];

			return result;
		};

		friend Matrix operator-(const Matrix &left, const double &right)// low efficient
		{
			printf("For debug, low efficient!\n");
			Matrix result(left.row, left.col);

			for (int i = 0; i < left.row; i++)
				for (int j = 0; j < left.col; j++)
					result.matrix[i + j * result.inc] = left.matrix[i + j * left.inc] - right;

			return result;
		};

		friend Matrix operator-(const double &left, const Matrix &right)// low efficient
		{
			printf("For debug, low efficient!\n");
			Matrix result(right.row, right.col);

			for (int i = 0; i < right.row; i++)
				for (int j = 0; j < right.col; j++)
					result.matrix[i + j * result.inc] = left - right.matrix[i + j * right.inc];

			return result;
		};

		friend Matrix operator*(const Matrix &left, const Matrix &right)// low efficient
		{
			printf("For debug, low efficient!\n");
			assert(left.col == right.row);
			Matrix result(left.row, right.col);
			for (int i = 0; i < right.col; i++)
				for (int j = 0; j < left.row; j++)
					for (int k = 0; k < left.col; k++)
						result.matrix[j + i * result.inc] += left.matrix[j + k * left.inc] * right.matrix[k + i * right.inc];

			return result;
		};

		friend Matrix operator*(const double &value, const Matrix &mat)// low efficient
		{
			printf("For debug, low efficient!\n");
			Matrix result(mat.row, mat.col);

			for (int i = 0; i < mat.row; i++)
				for (int j = 0; j < mat.col; j++)
					result.matrix[i + j * result.inc] = value * mat.matrix[i + j * mat.inc];

			return result;
		};

		friend Matrix operator*(const Matrix &mat, const double &value)// low efficient
		{
			printf("For debug, low efficient!\n");
			Matrix result(mat.row, mat.col);
			for (int i = 0; i < mat.row; i++)
				for (int j = 0; j < mat.col; j++)
					result.matrix[i + j * result.inc] = value * mat.matrix[i + j * mat.inc];

			return result;
		};

		friend Matrix operator/(const Matrix &mat, const double &value) // low efficient
		{
			printf("For debug, low efficient!\n");
			Matrix result(mat.row, mat.col);

			for (int i = 0; i < mat.row; i++)
				for (int j = 0; j < mat.col; j++)
					result.matrix[i + j * result.inc] = mat.matrix[i + j * mat.inc] / value;

			return result;
		};
		// end above low efficient functions ==========================================================

	public:
		/*Sparse square matrix times a dense matrix
		The sparse m by n matrix B is stored as three arrays:
			B for entries
			ir and jc use the same format as the index of mxgetir and mxgetjc in Matlab
			nzmax for the number of nonzero entries
		X is a n by p dense matrix */
		static void SPBtimesX(const double *B, const unsigned long long *ir, const unsigned long long *jc, integer nzmax, const double *X, integer m, integer n, integer p, double *result);
		
		/*a dense tall thin matrix times a sparse square matrix
		The sparse m by n matrix B is stored as three arrays:
		B for entries
		ir and jc use the same format as the index of mxgetir and mxgetjc in Matlab
		nzmax for the number of nonzero entries
		X is a p by m dense matrix */
		static void XtimesSPB(const double *X, const double *B, const unsigned long long *ir, const unsigned long long *jc, integer nzmax, integer p, integer m, integer n, double *result);

		// matrix multiplication
		static void DGEMM(double alpha, Matrix &M1, bool trans1, Matrix &M2, bool trans2, double beta, Matrix &result);
		static void CGEMM(doublecomplex alpha, Matrix &M1, bool trans1, Matrix &M2, bool trans2, doublecomplex beta, Matrix &result);

		/*Use dsyevx to compute the eigenvalue decomposition for a real symmetric matrix. This function
		computes all the eigenvalues and eigenvectors. If one needs only eigenvalues, please use dsyevx
		directly for efficient purpose. Note that the refered entries in S matrix will be destroyed after
		calling this function.*/
		static void EigenSymmetricM(char *UorL, Matrix &S, Matrix &eigenvalues, Matrix &eigenvectors);

		/* Compute matrix exponential exp(M) for a symmetric matrix M.
		Use EigenSymmetric function to compute the eigenvalue decomposition M = U D U^T.
		The matrix exponential is U exp(D) U^T. Note that the refered entries in S matrix
		will be destroyed after calling this function. If S and ExpS are the same variable,
		then the solution can be computed correctly and stored in ExpS.*/
		static void ExpSymmetricM(char *UorL, Matrix &S, Matrix &ExpS);

		/* Compute matrix exponential log(M) for a symmetric positive definite matrix M.
		Use EigenSymmetric function to compute the eigenvalue decomposition M = U D U^T.
		The matrix exponential is U log(D) U^T. Note that the refered entries in S matrix
		will be destroyed after calling this function. If S and ExpS are the same variable,
		then the solution can be computed correctly and stored in ExpS.*/
		static void LogSymmetricM(char *UorL, Matrix &S, Matrix &LogS);

		// solve the complex Sylevster equation A X + X B = C
		// On output: A and B are replaced by schur and negative schur forms and C is replace by the solution X.
		// A and B can be a same variable
		static void CSYL(Matrix &A, Matrix &B, Matrix &C);

		const Matrix &operator=(const Matrix &);

		double *matrix;
		integer row, col;
		integer inc;

		Matrix(const double *mat, integer r, integer c, integer i = 0)
		{
			row = r;
			col = c;
			inc = (i == 0) ? r : i;
			matrix = const_cast<double *> (mat);
			isfree = false;
		};

		~Matrix()
		{
			if (row == 0 || col == 0)
				return;
			if (isfree)
				delete[] matrix;
		};

	private:
		bool isfree;


		// below functions are used in the low efficient friend functions========================
		Matrix(const int r, const int c)
		{
			row = r;
			col = c;
			inc = r;
			matrix = new double[row];
			for (int i = 0; i < row * col; i++)
				matrix[i] = 0;
			isfree = true;
		};

		Matrix(Matrix &M)
		{
			row = M.row;
			col = M.col;
			matrix = new double[row];
			for (int i = 0; i < row * col; i++)
				matrix[i] = M.matrix[i];
		};
	};
}; /*end of ROPTLIB namespace*/
#endif
