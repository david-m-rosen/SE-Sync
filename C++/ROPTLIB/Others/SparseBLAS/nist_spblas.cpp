/*
This is the C++ version of Sparse BLAS downloaded from
http://math.nist.gov/spblas/
I slightly modified it such that it is compatible with
ROPTLIB.


--Wen Huang
*/

/*
*
* Sparse BLAS (Basic Linear Algebra Subprograms) Library
*
* A C++ implementation of the routines specified by the ANSI C 
* interface specification of the Sparse BLAS in the BLAS Technical 
* Forum Standard[1].   For details, see [2].
*
* Mathematical and Computational Sciences Division
* National Institute of Technology,
* Gaithersburg, MD USA
*
*
* [1] BLAS Technical Forum: www.netlib.org/blas/blast-forum/
* [2] I. S. Duff, M. A. Heroux, R. Pozo, "An Overview of the Sparse Basic
*     Linear Algebra Subprograms: The new standard of the BLAS Techincal
*     Forum,"  Vol. 28, No. 2, pp. 239-267,ACM Transactions on Mathematical 
*     Software (TOMS), 2002.
*
*
* DISCLAIMER:
*
* This software was developed at the National Institute of Standards and
* Technology (NIST) by employees of the Federal Government in the course
* of their official duties. Pursuant to title 17 Section 105 of the
* United States Code, this software is not subject to copyright protection
* and is in the public domain. NIST assumes no responsibility whatsoever for
* its use by other parties, and makes no guarantees, expressed or implied,
* about its quality, reliability, or any other characteristic.
*
*
*/


/* numeric is for accumulate() below */
#include <iostream>
#include <complex>
#include <numeric>
#include <vector>
#include <utility>
  /* pair defined here */

#include "Others/SparseBLAS/blas_sparse.h"




#ifdef SPBLAS_ERROR_FATAL
#include <cassert>
#define ASSERT_RETURN(x, ret_val) assert(x)
#define ERROR_RETURN(ret_val)  assert(0)
#else
#define ASSERT_RETURN(x, ret_val) {if (!(x)) return ret_val;}
#define ERROR_RETURN(ret_val) return ret_val
#endif


namespace NIST_SPBLAS
{

    
    /* dummy routines for real version of usdot to compile. */
    
    inline const double& conj(const double &x)
    {
        return x;
    }
    
    inline const float& conj(const float &x)
    { 
        return x;
    }
    

/**
   Generic sparse matrix (base) class: defines only the structure 
   (size, symmetry, etc.) and maintains state during construction, 
   but does not specify the actual nonzero values, or their type. 

*/
class Sp_mat
{
  private:

    int num_rows_;
    int num_cols_;
    int num_nonzeros_;

    /* ... */

    int void_;
    int nnew_;      /* avoid using "new" since it is a C++ keyword */
    int open_;
    int valid_;

    int unit_diag_ ;
    int complex_;
    int real_;
    int single_precision_;
    int double_precision_;
    int upper_triangular_;
    int lower_triangular_;
    int upper_symmetric_;
    int lower_symmetric_;
    int upper_hermitian_;
    int lower_hermitian_;
    int general_;

    int one_base_;


        /* optional block information */

    int Mb_;                /* matrix is partitioned into Mb x Nb blocks    */
    int Nb_;                /* otherwise 0, if regular (non-blocked) matrix */
    int k_;                 /* for constant blocks, each block is k x l     */
    int l_;                 /* otherwise 0, if variable blocks are used.   */

    int rowmajor_;          /* 1,if block storage is rowm major.  */
    int colmajor_;          /* 1,if block storage is column major. */

    /* unused optimization paramters */

    int opt_regular_;
    int opt_irregular_;
    int opt_block_;
    int opt_unassembled_;

    std::vector<int> K_; /* these are GLOBAL index of starting point of block     */
    std::vector<int> L_; /* i.e. block(i,j) starts at global location (K[i],L[i]) */
                    /* and of size (K[i+1]-K[i] x L[i+1]-L[i])               */

  public:

    Sp_mat(int M, int N) : 
      num_rows_(M),         /* default construction */
      num_cols_(N),
      num_nonzeros_(0),

      void_(0),
      nnew_(1),
      open_(0),
      valid_(0),

      unit_diag_(0),
      complex_(0),
      real_(0),
      single_precision_(0),
      double_precision_(0),
      upper_triangular_(0),
      lower_triangular_(0),
      upper_symmetric_(0),
      lower_symmetric_(0),
      upper_hermitian_(0),
      lower_hermitian_(0),
      general_(0),
      one_base_(0),
      Mb_(0),
      Nb_(0),
      k_(0),
      l_(0),
      rowmajor_(0),
      colmajor_(0),
      opt_regular_(0),
      opt_irregular_(1),
      opt_block_(0),
      opt_unassembled_(0),
      K_(),
      L_()
      {}


    int& num_rows()           { return num_rows_; }
    int& num_cols()           { return num_cols_; }
    int& num_nonzeros()         { return num_nonzeros_;}

    int num_rows() const        { return num_rows_; }
    int num_cols() const        { return num_cols_; }
    int num_nonzeros() const      { return num_nonzeros_;}

    int is_one_base() const     { return (one_base_ ? 1 : 0); }
    int is_zero_base() const    { return (one_base_ ? 0 : 1); }
    int is_void() const         { return void_; }
    int is_new() const          { return nnew_; }
    int is_open() const         { return open_; }
    int is_valid() const        { return valid_; }

    int is_unit_diag() const    { return unit_diag_; }
    int is_complex() const        { return complex_;}
    int is_real() const         { return real_;}
    int is_single_precision() const   { return single_precision_;}
    int is_double_precision() const   { return double_precision_;}
    int is_upper_triangular() const   { return upper_triangular_;}
    int is_lower_triangular() const   { return lower_triangular_;}
    int is_triangular() const     { return upper_triangular_ ||
                           lower_triangular_; }


    int is_lower_symmetric() const    { return lower_symmetric_; }
    int is_upper_symmetric() const    { return upper_symmetric_; }
    int is_symmetric() const      { return upper_symmetric_ ||
                           lower_symmetric_; }

    int is_lower_hermitian() const    { return lower_hermitian_; }
    int is_upper_hermitian() const    { return upper_hermitian_; }
    int is_hermitian() const  { return lower_hermitian_ || 
                                       upper_hermitian_; }
    int is_general() const { return !( is_hermitian() || is_symmetric()) ; }

    int is_lower_storage() const { return is_lower_triangular() ||
                                          is_lower_symmetric()  ||
                                          is_lower_hermitian() ; }

    int is_upper_storage() const { return is_upper_triangular() ||
                                          is_upper_symmetric()  ||
                                          is_upper_hermitian() ; }


    int is_opt_regular() const { return opt_regular_; }
    int is_opt_irregular() const { return opt_irregular_; }
    int is_opt_block() const { return opt_block_;} 
    int is_opt_unassembled() const { return opt_unassembled_;}

    int K(int i) const { return (k_ ? i*k_ : K_[i] ); }
    int L(int i) const { return (l_ ? i*l_ : L_[i] ); }

    int is_rowmajor() const { return rowmajor_; }
    int is_colmajor() const { return colmajor_; }


    void set_one_base()   { one_base_ = 1; }
    void set_zero_base()  { one_base_ = 0; }


    void set_void()       { void_ = 1;  nnew_ = open_ =  valid_ = 0;}
    void set_new()        { nnew_ = 1;  void_ = open_ =  valid_ = 0;}
    void set_open()       { open_ = 1;  void_ = nnew_  = valid_ = 0;}
    void set_valid()      { valid_ = 1; void_ = nnew_ =  open_ = 0; }


    void set_unit_diag()    { unit_diag_ = 1;}
    void set_complex()        {complex_ = 1; }
    void set_real()         { real_ = 1; }
    void set_single_precision()   { single_precision_ = 1; }
    void set_double_precision()   { double_precision_ = 1; }
    void set_upper_triangular()   { upper_triangular_ = 1; }
    void set_lower_triangular()   { lower_triangular_ = 1; }
    void set_upper_symmetric()  { upper_symmetric_ = 1; }
    void set_lower_symmetric()  { lower_symmetric_ = 1; }
    void set_upper_hermitian()  { upper_hermitian_ = 1; }
    void set_lower_hermitian()  { lower_hermitian_ = 1; }

  
    void set_const_block_parameters(int Mb, int Nb, int k, int l)
    {
      Mb_ = Mb;
      Nb_ = Nb;
      k_ = k;
      l_ = l;
    }


    void set_var_block_parameters(int Mb, int Nb, const int *k, const int *l)
    {
      Mb_ = Mb;
      Nb_ = Nb;
      k_ = 0;
      l_ = 0;

      K_.resize(Mb+1);
      K_[0] = 0;
      for (int i=0; i<Mb; i++)
        K_[i+1] = k[i] + K_[i];
      
      L_.resize(Nb+1);
      L_[0] = 0;
      for (int j=0; j<Mb; j++)
        K_[j+1] = k[j] + K_[j];
    }


    virtual int end_construction()
    {
      if (is_open() || is_new())
      {
        set_valid();

        return 0;
      }
      else
        ERROR_RETURN(1);
    }
    virtual void print() const;

    virtual void destroy() {};

    virtual ~Sp_mat() {};

};


template <class T>
class TSp_mat : public Sp_mat
{
  private:
    std::vector< std::vector< std::pair<T, int> > > S;
    std::vector<T> diag;                 /* optional diag if matrix is
                        triangular. Created
                        at end_construction() phase */
  private:


    inline T sp_dot_product( const std::vector< std::pair<T, int> > &r, 
        const T* x, int incx ) const
    {
        T sum(0);

        if (incx == 1)
        {
          for ( typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
            p < r.end(); p++)
          {
            //sum = sum + p->first * x[p->second];
            sum += p->first * x[p->second];
          }
        }
        else /* incx != 1 */
        {
          for ( typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
            p < r.end(); p++)
           {
            //sum = sum + p->first * x[p->second * incx];
            sum += p->first * x[p->second * incx];
            }
        }
    
        return sum;
  }

    inline T sp_conj_dot_product( const std::vector< std::pair<T, int> > &r, 
        const T* x, int incx ) const
    {
        T sum(0);

        if (incx == 1)
        {
          for ( typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
            p < r.end(); p++)
          {
            sum += conj(p->first) * x[p->second];
          }
        }
        else /* incx != 1 */
        {
          for ( typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
            p < r.end(); p++)
           {
            //sum = sum + p->first * x[p->second * incx];
            sum += conj(p->first) * x[p->second * incx];
            }
        }
    
        return sum;
  }


  inline void sp_axpy( const T& alpha, const std::vector< std::pair<T,int> > &r, 
      T*  y, int incy) const
  {
    if (incy == 1)
    {
      for (typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
          p < r.end(); p++)
       y[p->second] += alpha * p->first;  

    }
    else /* incy != 1 */
    {
    for (typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
        p < r.end(); p++)
      y[incy * p->second] += alpha * p->first;  
    } 
  } 

  inline void sp_conj_axpy( const T& alpha, const std::vector< std::pair<T,int> > &r, 
      T*  y, int incy) const
  {
    if (incy == 1)
    {
      for (typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
          p < r.end(); p++)
       y[p->second] += alpha * conj(p->first);  

    }
    else /* incy != 1 */
    {
    for (typename std::vector< std::pair<T,int> >::const_iterator p = r.begin(); 
        p < r.end(); p++)
      y[incy * p->second] += alpha * conj(p->first);  
    } 
  } 

  void mult_diag(const T& alpha, const T* x, int incx, T* y, int incy) 
      const
  {
    const T* X = x;
    T* Y = y;
    typename std::vector<T>::const_iterator d= diag.begin();
    for ( ; d < diag.end(); X+=incx, d++, Y+=incy)
    {
      *Y += alpha * *d * *X;
    }
  }

  void mult_conj_diag(const T& alpha, const T* x, int incx, T* y, int incy) 
      const
  {
    const T* X = x;
    T* Y = y;
    typename std::vector<T>::const_iterator d= diag.begin();
    for ( ; d < diag.end(); X+=incx, d++, Y+=incy)
    {
      *Y += alpha * conj(*d) * *X;
    }
  }


  void nondiag_mult_vec(const T& alpha, const T* x, int incx, 
      T* y, int incy) const
  {

    int M = num_rows();

    if (incy == 1)
    {
      for (int i=0; i<M; i++)
        y[i] += alpha * sp_dot_product(S[i], x, incx);
    }
    else
    {
      for (int i=0; i<M; i++)
        y[i * incy] += alpha * sp_dot_product(S[i], x, incx);
    }
  }

  void nondiag_mult_vec_conj(const T& alpha, const T* x, int incx, 
      T* y, int incy) const
  {

    int M = num_rows();

    if (incy == 1)
    {
      for (int i=0; i<M; i++)
        y[i] += alpha * sp_conj_dot_product(S[i], x, incx);
    }
    else
    {
      for (int i=0; i<M; i++)
        y[i * incy] += alpha * sp_conj_dot_product(S[i], x, incx);
    }
  }

  void nondiag_mult_vec_transpose(const T& alpha, const T* x, int incx, 
      T* y, int incy) const
  {
    /* saxpy: y += (alpha * x[i]) row[i]  */

    int M = num_rows();
    const T* X = x;
    for (int i=0; i<M; i++, X += incx)
      sp_axpy( alpha * *X, S[i], y, incy);
  }

  void nondiag_mult_vec_conj_transpose(const T& alpha, const T* x, int incx, 
      T* y, int incy) const
  {
    /* saxpy: y += (alpha * x[i]) row[i]  */

    int M = num_rows();
    const T* X = x;
    for (int i=0; i<M; i++, X += incx)
      sp_conj_axpy( alpha * *X, S[i], y, incy);
  }

  void mult_vec(const T& alpha, const T* x, int incx, T* y, int incy) 
      const
  {
    nondiag_mult_vec(alpha, x, incx, y, incy);

    if (is_triangular() || is_symmetric())
      mult_diag(alpha, x, incx, y, incy);

    if (is_symmetric())
      nondiag_mult_vec_transpose(alpha, x, incx, y, incy);
  }


  void mult_vec_transpose(const T& alpha, const T* x, int incx, T* y, 
      int incy) const
  {

    nondiag_mult_vec_transpose(alpha, x, incx, y, incy);

    if (is_triangular() || is_symmetric())
      mult_diag(alpha, x, incx, y, incy);

    if (is_symmetric())
      nondiag_mult_vec(alpha, x, incx, y, incy);
  }


  void mult_vec_conj_transpose(const T& alpha, const T* x, int incx, T* y, 
      int incy) const
  {

    nondiag_mult_vec_conj_transpose(alpha, x, incx, y, incy);

    if (is_triangular() || is_symmetric())
      mult_conj_diag(alpha, x, incx, y, incy);

    if (is_symmetric())
      nondiag_mult_vec_conj(alpha, x, incx, y, incy);
  }






      
  int triangular_solve(T alpha, T* x, int incx ) const
  {
    if (alpha == (T) 0.0)
      ERROR_RETURN(1);

    if ( ! is_triangular() )
      ERROR_RETURN(1);

    int N = num_rows();

    if (is_lower_triangular())
    {
        for (int i=0, ii=0; i<N; i++, ii += incx)
        {
            x[ii] = (x[ii] - sp_dot_product(S[i], x, incx)) / diag[i];
        }
       if (alpha != (T) 1.0)
       {
        for (int i=0, ii=0; i<N; i++, ii += incx)
            x[ii] /= alpha; 
       }
    }
    else if (is_upper_triangular())
    {

      for (int i=N-1, ii=(N-1)*incx ;   0<=i ;    i--, ii-=incx)
      {
         x[ii] = (x[ii] - sp_dot_product(S[i],x, incx)) / diag[i];
      }
      if (alpha != (T) 1.0)
      {
        for (int i=N-1, ii=(N-1)*incx ;   0<=i ;    i--, ii-=incx)
          x[ii] /= alpha; 
      }

    }
    else
        ERROR_RETURN(1);

    return 0;
  }

  int transpose_triangular_solve(T alpha, T* x, int incx) const
  {
    if ( ! is_triangular())
      return -1;

    int N = num_rows();

    if (is_lower_triangular())
    {

      for (int j=N-1, jj=(N-1)*incx; 0<=j; j--, jj -= incx)
      {
        x[jj] /= diag[j] ;
        sp_axpy( -x[jj], S[j], x, incx);
      }
      if (alpha != (T) 1.0)
      {
        for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
          x[jj] /= alpha;
      }
    }
    else if (is_upper_triangular())
    {
      
      for (int j=0, jj=0; j<N; j++, jj += incx)
      {
        x[jj] /= diag[j];
        sp_axpy(- x[jj], S[j], x, incx);
      }
      if (alpha != (T) 1.0)
      {
        for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
          x[jj] /= alpha;
      }
    }
    else
         ERROR_RETURN(1);

    return 0;
  }



  int transpose_triangular_conj_solve(T alpha, T* x, int incx) const
  {
    if ( ! is_triangular())
      return -1;

    int N = num_rows();

    if (is_lower_triangular())
    {

      for (int j=N-1, jj=(N-1)*incx; 0<=j; j--, jj -= incx)
      {
        x[jj] /= conj(diag[j]) ;
        sp_conj_axpy( -x[jj], S[j], x, incx);
      }
      if (alpha != (T) 1.0)
      {
        for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
          x[jj] /= alpha;
      }
    }
    else if (is_upper_triangular())
    {
      
      for (int j=0, jj=0; j<N; j++, jj += incx)
      {
        x[jj] /= conj(diag[j]);
        sp_conj_axpy(- x[jj], S[j], x, incx);
      }
      if (alpha != (T) 1.0)
      {
        for (int jj=(N-1)*incx; 0<=jj; jj -=incx)
          x[jj] /= alpha;
      }
    }
    else
         ERROR_RETURN(1);

    return 0;
  }




 public:

  inline T& val(std::pair<T, int> &VP) { return VP.first; }
  inline int& col_index(std::pair<T,int> &VP) { return VP.second; } 

  inline const T& val(std::pair<T, int> const &VP) const { return VP.first; }
  inline int col_index(std::pair<T,int> const &VP) const { return VP.second; } 


  TSp_mat( int M, int N) : Sp_mat(M,N), S(M), diag() {}

  void destroy()
  {
    // set std::vector sizes to zero
    (std::vector<T>(0)).swap(diag);
    (std::vector< std::vector< std::pair<T, int> > > (0) ).swap(S);
  }

/**

    This function is the entry point for all of the insert routines in 
    this implementation.  It fills the sparse matrix, one entry at a time.
    If matrix is declared unit_diagonal, then inserting any diagonal
    values is ignored.  If it is symmetric (upper/lower) or triangular
    (upper/lower) inconsistent values are not caught.  (That is, entries
    into the upper region of a lower triangular matrix is not reported.)

    [NOTE: the base is determined at the creation phase, and can be determined
    by testing whether  BLAS_usgp(A, blas_one_base) returns 1.  If it returns 0,
    then offsets are zero based.]

    @param val  the numeric value of entry A(i,j)
    @param i  the row index of A(i,j)  
    @param j  the column index of A(i,j)

    @return 0 if succesful, 1 otherwise
*/  
  int insert_entry(T val, int i, int j)
  {
    if (is_one_base())        
    {
      i--;
      j--;
    }

    /* make sure the indices are in range */
    ASSERT_RETURN(i >= 0, 1);
    ASSERT_RETURN(i < num_rows(), 1);
    ASSERT_RETURN(j >= 0, 1);
    ASSERT_RETURN(j < num_cols(), 1);

    /* allocate space for the diagonal, if this is the first time
     * trying to insert values.
    */
    if (is_new())
    {
      set_open();
      
      if (is_triangular() || is_symmetric())
      {
        diag.resize(num_rows());

        if (is_unit_diag())
        {
          for (unsigned int ii=0; ii< diag.size(); ii++)
              diag[ii] = T(1.0); 
        }
        else
        {
          for (unsigned int ii=0; ii< diag.size(); ii++)
              diag[ii] = (T) 0.0; 
        }
      }

    }
    if (is_open())
    {

      if (i==j && (is_triangular() || is_symmetric() || is_hermitian()) )
      {
        if (!is_unit_diag())
        {
          diag[i] += val;
        }
        else /* if unit diagonal */
        {
          if (val != (T) 1) 
            ERROR_RETURN(0);    /* tries to insert non-unit diagonal */
        }

        if (is_upper_storage() && i > j)
            ERROR_RETURN(0);    /* tries to fill lower-triangular region */
        else 
        
          if (is_lower_storage() && i < j)
            ERROR_RETURN(0);  /* tries to fill upper-triangular region */

      }
      else
      {
          S[i].push_back( std::make_pair(val, j) );
      }

      num_nonzeros() ++;
    }


    return 0;
  }

  int insert_entries( int nz, const T* Val, const int *I, const int *J)
  {
    for (int i=0; i<nz; i++)
    {
      insert_entry(Val[i], I[i], J[i]) ;
    }
    return 0;

  }

  int insert_row(int k, int nz, const T* Val, const int *J)
  {
    for (int i=0; i<nz; i++)
      insert_entry(Val[i], k, J[i]);  
    return 0;
  }

  int insert_col(int k, int nz, const T* Val, const int *I)
  {
    for (int i=0; i<nz; i++)
      insert_entry(Val[i], I[i], k);  
    return 0;
  }

  int insert_block(const T* Val, int row_stride, 
        int col_stride, int bi, int bj)
  {
    /* translate from block index to global indices */
    int Iend = K(bi+1);
    int Jend = L(bj+1);
    for (int i=K(bi), r=0; i<Iend; i++, r += row_stride)
      for (int j=L(bi); j<Jend; j++, r += col_stride)
        insert_entry( Val[r], i, j );

    return 0;
  }

  int end_construction()
  {
    return Sp_mat::end_construction();
  }



  int usmv(enum blas_trans_type transa, const T& alpha, const  T* x , int incx, 
    T* y, int incy) const
  {
    
  ASSERT_RETURN(is_valid(), -1);

  if (transa == blas_no_trans)
    mult_vec(alpha, x, incx, y, incy);
  else
  if (transa == blas_conj_trans)
    mult_vec_conj_transpose(alpha, x, incx, y, incy);
  else
  if ( transa == blas_trans)
    mult_vec_transpose(alpha, x, incx, y, incy);
  else
    ERROR_RETURN(1);
    
  
    return 0;

  }


  int usmm(enum blas_order_type ordera, enum blas_trans_type transa, 
    int nrhs, const T& alpha, const  T* b, int ldb, T* C, int ldC) const
  {
    if (ordera == blas_rowmajor)
    {
      /* for each column of C, perform a mat_vec */
      for (int i=0; i<nrhs; i++)
      {
        usmv( transa, alpha, &b[i], ldb, &C[i], ldC );
      }
      return 0;
    }
    else
    if (ordera == blas_colmajor)
    {
      /* for each column of C, perform a mat_vec */
      for (int i=0; i<nrhs; i++)
      {
        usmv( transa, alpha, &b[i*ldb], 1, &C[i*ldC], 1 );
      }
      return 0;
    }
    else
      ERROR_RETURN(1);
  }

  int ussv( enum blas_trans_type transa, const T& alpha,  T* x, int incx) const
  {
      if (transa == blas_trans)
        return transpose_triangular_solve(alpha, x, incx);
      else 
      if (transa == blas_conj_trans)
        return transpose_triangular_conj_solve(alpha, x, incx);
      else
      if (transa == blas_no_trans)
        return triangular_solve(alpha, x, incx);
      else
        ERROR_RETURN(1);


  }

  int ussm( enum blas_order_type ordera, enum blas_trans_type transa, int nrhs,
      const T& alpha, T* C, int ldC) const
  {
    if (ordera == blas_rowmajor)
    {
      /* for each column of C, perform a usmv */
      for (int i=0; i<nrhs; i++)
      {
        ussv( 
            transa, alpha, &C[i], ldC );
      }
      return 0;
    }
    else
    if (ordera == blas_colmajor)
    {
      /* for each column of C, perform a mat_vec */
      for (int i=0; i<nrhs; i++)
      {
        ussv( transa, alpha, &C[i*ldC], 1 );
      }
      return 0;
    }
    else
      ERROR_RETURN(1);
  } 



  void print() const
  {
    Sp_mat::print();  /* print matrix header info */

    /* if there is actual data, print out contents */
    for (int i=0; i<num_rows(); i++)
      for (unsigned int j=0; j< S[i].size(); j++)
        std::cout << i << "    " << col_index(S[i][j]) <<
              "        "  << val(S[i][j]) << "\n";

    /* if matrix is triangular, print out diagonals */
    if (is_upper_triangular() || is_lower_triangular())
    {
      for (unsigned int i=0; i< diag.size(); i++)
        std::cout << i << "    " << i << "     " << diag[i] << "\n";
    }
  }

    
};




typedef TSp_mat<float> FSp_mat;
typedef TSp_mat<double> DSp_mat;
typedef TSp_mat<std::complex<float> > CSp_mat;
typedef TSp_mat<std::complex<double> > ZSp_mat;



void table_print();
void print(int A);

}
/* namespace */


namespace NIST_SPBLAS
{
static std::vector<Sp_mat *> Table;
static unsigned int Table_active_matrices = 0;
int Table_insert(Sp_mat* S);
int Table_remove(unsigned int i);

/* 
  finds an empty slot in global sparse marix table, and fills it with
  given entry.  Returns -1 if no spot found, or table is corrupt.
*/
int Table_insert(Sp_mat* S)
{
  if (Table_active_matrices <= Table.size())
  {
    Table.push_back(S);
    Table_active_matrices++;
    return static_cast<int> (Table.size() - 1);
  }
  else
  {
    /* there is an available slot; find it. */
    for (unsigned int i=0; i<Table.size(); i++)
    {
      if (Table[i] == NULL)
      {
        Table[i] = S;
        Table_active_matrices++;
        return i;
      }
    }
  }

  return -1;
}

/* 
  removes an exisiting sparse matrix from global table.  Returns 0, if
  successfull, 1 otherwise.
*/
int Table_remove(unsigned int i)
{
  if (i < Table.size() && Table[i] != NULL)
  {
    Table[i] = NULL;
    Table_active_matrices--;
    return 0;
  }
  else
    return -1;
}



void Sp_mat::print() const
{
  

  std::cout << "State : " <<
    (is_void() ? "void" : 
     is_new()  ? "new" :  
     is_open() ? "open" :
     is_valid() ? "valid" : "unknown") << "\n";

  std::cout << "M = " <<  num_rows() <<  "  N = " << num_cols() <<
        "  nz = " << num_nonzeros() << "\n";

#define yesno(exp) ( (exp) ? "yes" : "no" )

  std::cout << "real: "     << yesno(is_real()) << "\n";
  std::cout << "complex: "  << yesno(is_complex()) << "\n";
  std::cout << "double "    << yesno(is_double_precision()) << "\n";
  std::cout << "single "    << yesno(is_single_precision()) << "\n";

  std::cout << "upper_triangular: " << yesno(is_upper_triangular()) << "\n";
  std::cout << "lower_triangular: " << yesno(is_lower_triangular()) << "\n";

  std::cout << "regular:    " << yesno(is_opt_regular()) << "\n";
  std::cout << "irregular:  " << yesno(is_opt_irregular()) << "\n";
  std::cout << "block:      " << yesno(is_opt_block()) << "\n";
  std::cout << "unassembled:" << yesno(is_opt_unassembled()) << "\n";
      
#undef yesno
}

void table_print()
{
  std::cout << "Table has " << Table.size() << " element(s). \n";
  for (unsigned int i=0; i< Table.size(); i++)
  {
    if (Table[i] != 0)
    {
      std::cout << "***** Table[" << i << "]: \n";
      Table[i]->print();
      std::cout << "\n\n";
    }
  }
}

void print(int A)
{
  std::cout << "\n";
  Table[A]->print();
  std::cout << "\n";
}



}
/* namespace NIST_SPBLAS */

using namespace NIST_SPBLAS;





/* Level 1 */

/* these macros are useful for creating some consistency between the 
   various precisions and floating point types.
*/

typedef float    FLOAT;
typedef double   DOUBLE;
typedef std::complex<float> COMPLEX_FLOAT;
typedef std::complex<double> COMPLEX_DOUBLE;



typedef float          SPBLAS_FLOAT_IN;
typedef double         SPBLAS_DOUBLE_IN;
typedef const void *   SPBLAS_COMPLEX_FLOAT_IN;
typedef const void *   SPBLAS_COMPLEX_DOUBLE_IN;

typedef float *  SPBLAS_FLOAT_OUT;
typedef double * SPBLAS_DOUBLE_OUT;
typedef void *   SPBLAS_COMPLEX_FLOAT_OUT;
typedef void *   SPBLAS_COMPLEX_DOUBLE_OUT;



typedef float *  SPBLAS_FLOAT_IN_OUT;
typedef double * SPBLAS_DOUBLE_IN_OUT;
typedef void *   SPBLAS_COMPLEX_FLOAT_IN_OUT;
typedef void *   SPBLAS_COMPLEX_DOUBLE_IN_OUT;

typedef const float *  SPBLAS_VECTOR_FLOAT_IN;
typedef const double * SPBLAS_VECTOR_DOUBLE_IN;
typedef const void *   SPBLAS_VECTOR_COMPLEX_FLOAT_IN;
typedef const void *   SPBLAS_VECTOR_COMPLEX_DOUBLE_IN;

typedef float *  SPBLAS_VECTOR_FLOAT_OUT;
typedef double * SPBLAS_VECTOR_DOUBLE_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_FLOAT_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_DOUBLE_OUT;

typedef float *  SPBLAS_VECTOR_FLOAT_IN_OUT;
typedef double * SPBLAS_VECTOR_DOUBLE_IN_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT;
typedef void *   SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT;



#define SPBLAS_TO_FLOAT_IN(x)   x
#define SPBLAS_TO_DOUBLE_IN(x)  x
#define SPBLAS_TO_COMPLEX_FLOAT_IN(x) \
        (* reinterpret_cast<const std::complex<float> *>(x))
#define SPBLAS_TO_COMPLEX_DOUBLE_IN(x)  \
        (* reinterpret_cast<const std::complex<double> *>(x))

#define SPBLAS_TO_FLOAT_OUT(x)  x
#define SPBLAS_TO_DOUBLE_OUT(x) x
#define SPBLAS_TO_COMPLEX_FLOAT_OUT(x)  reinterpret_cast<std::complex<float> *>(x)
#define SPBLAS_TO_COMPLEX_DOUBLE_OUT(x) reinterpret_cast<std::complex<double> *>(x)  

#define SPBLAS_TO_FLOAT_IN_OUT(x)   x
#define SPBLAS_TO_DOUBLE_IN_OUT(x)  x
#define SPBLAS_TO_COMPLEX_FLOAT_IN_OUT(x)  reinterpret_cast<std::complex<float> *>(x)
#define SPBLAS_TO_COMPLEX_DOUBLE_IN_OUT(x) reinterpret_cast<std::complex<double>*>(x)  


#define SPBLAS_TO_VECTOR_DOUBLE_IN(x)   x 
#define SPBLAS_TO_VECTOR_FLOAT_IN(x)  x 
#define SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN(x) \
                          reinterpret_cast<const std::complex<float>*>(x)
#define SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN(x) \
                          reinterpret_cast<const std::complex<double>*>(x)

#define SPBLAS_TO_VECTOR_DOUBLE_OUT(x)  x 
#define SPBLAS_TO_VECTOR_FLOAT_OUT(x)   x 
#define SPBLAS_TO_VECTOR_COMPLEX_FLOAT_OUT(x) \
                          reinterpret_cast<std::complex<float>*>(x)
#define SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_OUT(x) \
                          reinterpret_cast<std::complex<double>*>(x)

#define SPBLAS_TO_VECTOR_DOUBLE_IN_OUT(x)   x 
#define SPBLAS_TO_VECTOR_FLOAT_IN_OUT(x)  x 
#define SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT(x) \
                          reinterpret_cast<std::complex<float>*>(x)
#define SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT(x) \
                          reinterpret_cast<std::complex<double>*>(x)






#define BLAS_FLOAT_NAME(routine_name) BLAS_s##routine_name
#define BLAS_DOUBLE_NAME(routine_name) BLAS_d##routine_name
#define BLAS_COMPLEX_FLOAT_NAME(routine_name) BLAS_c##routine_name
#define BLAS_COMPLEX_DOUBLE_NAME(routine_name) BLAS_z##routine_name


#define TSp_MAT_SET_FLOAT(A) {A->set_single_precision(); A->set_real();}
#define TSp_MAT_SET_DOUBLE(A) {A->set_double_precision(); A->set_real();}
#define TSp_MAT_SET_COMPLEX_FLOAT(A) {A->set_single_precision(); A->set_complex();}
#define TSp_MAT_SET_COMPLEX_DOUBLE(A) {A->set_double_precision(); A->set_complex();}


        /*------------------------------------*/
        /* Non-precision Sparse BLAS routines */
        /*------------------------------------*/

/* -------- */
/*  USSP()  */
/* -------- */
int BLAS_ussp(blas_sparse_matrix A, int pname)
{
    Sp_mat *S = Table[A];


                                /* Note: these are returns, in the case */
                                /* statement, so "break" is not needed.  */
    switch (pname)
    {
        case (blas_zero_base) : S->set_zero_base(); break;
        case (blas_one_base)  : S->set_one_base(); break;

        case (blas_unit_diag) : S->set_unit_diag(); break;
        case (blas_complex)   : S->set_complex(); break;
        case (blas_real)      : S->set_real(); break;

        case (blas_single_precision) : S->set_single_precision(); break;
        case (blas_double_precision) : S->set_double_precision(); break;


        case (blas_lower_triangular) : S->set_lower_triangular(); break;
        case (blas_upper_triangular) : S->set_upper_triangular(); break;

        case (blas_lower_symmetric) : S->set_lower_symmetric(); break;
        case (blas_upper_symmetric) : S->set_upper_symmetric(); break;

        case (blas_lower_hermitian) : S->set_lower_hermitian(); break;
        case (blas_upper_hermitian) : S->set_upper_hermitian(); break;


                                        /* optimizations not used */
        case (blas_regular )        :
        case (blas_irregular)       : 
        case (blas_block)           :
        case (blas_unassembled)     :    return 0;

        default:                return -1;  /* invalid property */
    }

    return 0;

}

/* -------- */
/*  USGP()  */
/* -------- */

int BLAS_usgp(blas_sparse_matrix A, int pname)
{
  Sp_mat *S = Table[A];

    switch (pname)
    {
        case (blas_num_rows)          : return S->num_rows(); 
        case (blas_num_cols)          : return S->num_cols(); 
        case (blas_num_nonzeros)      : return S->num_nonzeros(); 

        case (blas_complex)           : return S->is_complex();
        case (blas_real)              : return S->is_real();
        case (blas_single_precision)  : return S->is_single_precision();
        case (blas_double_precision)  : return S->is_double_precision();


        case (blas_lower_triangular) : return S->is_lower_triangular(); break;
        case (blas_upper_triangular) : return S->is_upper_triangular(); break;

        case (blas_general)           : return S->is_general();
        case (blas_symmetric)         : return S->is_symmetric();
        case (blas_hermitian)         : return S->is_hermitian();

        case (blas_zero_base) : return S->is_zero_base(); 
        case (blas_one_base) : return  S->is_one_base(); 


        case (blas_rowmajor)          : return S->is_rowmajor();
        case (blas_colmajor)          : return S->is_colmajor();
        case (blas_new_handle)        : return S->is_new();
        case (blas_valid_handle)      : return S->is_valid();
        case (blas_open_handle)       : return S->is_open();
        case (blas_invalid_handle)    : return S->is_void();

        case (blas_regular)           : return S->is_opt_regular();
        case (blas_irregular)         : return S->is_opt_irregular();
        case (blas_block)             : return S->is_opt_block();
        case (blas_unassembled)       : return S->is_opt_unassembled();

        default:                return -1;  /* invalid property */
    }
} 

/* -------- */
/*  USDS()  */
/* -------- */
int BLAS_usds(int A)
{
  Sp_mat *S  = Table[A];

  S->destroy();
  Table_remove(A);

  return 0;
}




/* --------------------------- */
/*  Level 1 generic routines   */
/* --------------------------- */

template <class T>
void BLAS_xusdot( enum blas_conj_type conj_flag, int nz,
    const T *x,  const int *index,  const T *y, int incy,
    T *r, enum blas_base_type index_base)
{

  T t(0);

  if (index_base == blas_one_base)
    y -= incy;

  if (conj_flag == blas_no_conj)
  {
    for (int i=0; i<nz; i++)
      t += x[i] * y[index[i]*incy];
  }
  else
    for (int i=0; i<nz; i++)
      t += conj(x[i]) * y[index[i]*incy];


  *r = t;
}


template <class T>
void BLAS_xusaxpy(int nz, T alpha, const T *x,
    const int *index, T *y, int incy,
    enum blas_base_type index_base)
{

  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
  {
//     y[index[i]*incy] +=  (alpha * x[i]);

  }
}

template <class T>
void BLAS_xusga( int nz, const T *y, int incy, T *x, const int *indx,
              enum blas_base_type index_base )
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
    x[i] = y[indx[i]*incy];

}

template <class T>
void BLAS_xusgz( int nz, T *y, int incy, T *x, const int *indx,
              enum blas_base_type index_base )
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
  {
    x[i] = y[indx[i]*incy];
    y[indx[i]*incy] = (T) 0.0;
  }

}


template <class T>
void BLAS_xussc(int nz, const T *x, T *y, int incy, const int *index,
    enum blas_base_type index_base)
{
  if (index_base == blas_one_base)
    y -= incy;

  for (int i=0; i<nz; i++)
    y[index[i]*incy] = x[i];

}

/* --------------------------- */
/* Level 2&3 generic precision */
/* --------------------------- */


template <class T>
int BLAS_xuscr_insert_entry(blas_sparse_matrix A,  const T& val, int i, int j)
{
  return ((TSp_mat<T> *)Table[A])->insert_entry(val, i, j);

}

template <class T> 
int BLAS_xuscr_insert_entries(blas_sparse_matrix A, int nz, const T* Val, 
      const int* I, const int *J)
{
  return ((TSp_mat<T>*) Table[A])->insert_entries(nz, Val, I, J);
}


template <class T> 
int BLAS_xuscr_insert_col(blas_sparse_matrix A, int j, int nz, const T* Val, 
      const int* indx)
{
  return ((TSp_mat<T>*) Table[A])->insert_col(j, nz, Val, indx);
}

template <class T> 
int BLAS_xuscr_insert_row(blas_sparse_matrix A, int i, int nz, const T* Val, 
      const int* indx)
{
  return ((TSp_mat<T>*) Table[A])->insert_row(i, nz, Val, indx);
}

template <class T> 
int BLAS_xuscr_insert_clique(blas_sparse_matrix A, int k, int l, const T* Val, 
      const int row_stride, const int col_stride, const int *indx, 
      const int *jndx)
{
  return ((TSp_mat<T>*) Table[A])->insert_clique(k, l, Val, row_stride, 
            col_stride, indx, jndx);
}

template <class T> 
int BLAS_xuscr_insert_block(blas_sparse_matrix A, const T* Val, 
      const int row_stride, const int col_stride, int bi, int bj )
{
  return ((TSp_mat<T>*) Table[A])->insert_block(Val, 
        row_stride, col_stride, bi, bj); 
}

inline int BLAS_xuscr_end(blas_sparse_matrix A)
{
  return (Table[A])->end_construction();
}


template <class T>
int BLAS_xusmv(enum blas_trans_type transa, const T& alpha, 
    blas_sparse_matrix A, const T *x, int incx, T *y, int incy )
{
  TSp_mat<T> *M = (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->usmv(transa, alpha, x, incx, y, incy);
}


template <class T>
int BLAS_xusmm(enum blas_order_type ordera, enum blas_trans_type transa, 
    int nrhs, const T& alpha, blas_sparse_matrix A, 
    const T *B, int ldB, T* C, int ldC)
{
  TSp_mat<T> *M = (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->usmm(ordera, transa, nrhs, alpha, B, ldB, C, ldC);
}


template <class T>
int BLAS_xussv(enum blas_trans_type transa, const T& alpha, 
    blas_sparse_matrix A, T *x, int incx)
{
  TSp_mat<T> *M = 
        (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->ussv(transa, alpha, x, incx);
}

template <class T>
int BLAS_xussm(enum blas_order_type orderA, enum blas_trans_type transa, 
    int nrhs, const T& alpha, blas_sparse_matrix A, T *C, int ldC)
{
  TSp_mat<T> *M = 
        (TSp_mat<T> *) Table[A];

  ASSERT_RETURN(M->is_valid(), 1);

  return M->ussm(orderA, transa, nrhs, alpha, C, ldC);
}

/*   --- end of generic rouintes ---- */


/*********/
/*  ----  double  Level 1 rouintes ----- */
/*********/

 void BLAS_DOUBLE_NAME(usdot)( 
    enum blas_conj_type conj_flag, 
    int  nz, 
    SPBLAS_VECTOR_DOUBLE_IN x,  
    const int *index,  
    SPBLAS_VECTOR_DOUBLE_IN y, 
    int incy, 
    SPBLAS_DOUBLE_OUT r, 
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz, 
          SPBLAS_TO_VECTOR_DOUBLE_IN( x ), index, 
          SPBLAS_TO_VECTOR_DOUBLE_IN( y ), incy, 
          SPBLAS_TO_DOUBLE_OUT( r ), index_base);
}

 void BLAS_DOUBLE_NAME(usaxpy)(
    int nz, 
    SPBLAS_DOUBLE_IN alpha, 
    SPBLAS_VECTOR_DOUBLE_IN x, 
    const int *index, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT y, 
    int incy, 
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_DOUBLE_IN( alpha ), 
      SPBLAS_TO_VECTOR_DOUBLE_IN( x ), index, 
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( y ), 
      incy, index_base);
}

 void BLAS_DOUBLE_NAME(usga)( 
    int nz, 
    SPBLAS_VECTOR_DOUBLE_IN y, 
    int incy, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_DOUBLE_IN( y ), incy, 
        SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( x ), indx, index_base);

}

 void BLAS_DOUBLE_NAME(usgz)( 
    int nz, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT y, 
    int incy, 
    SPBLAS_VECTOR_DOUBLE_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_DOUBLE_IN_OUT(y ), incy, SPBLAS_TO_DOUBLE_OUT(  x ), 
    indx, index_base);
}


 void BLAS_DOUBLE_NAME(ussc)(
    int nz, 
    SPBLAS_VECTOR_DOUBLE_IN x, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT y, 
    int incy, 
    const int *index, 
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_DOUBLE_IN( x ), 
  SPBLAS_TO_DOUBLE_IN_OUT( y ), incy, index, 
  index_base);
}

/*  DOUBLE Level 2/3 creation routines */

int BLAS_DOUBLE_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<DOUBLE> *A = new TSp_mat<DOUBLE>(M, N);
  TSp_MAT_SET_DOUBLE(A);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_DOUBLE_NAME(uscr_block_begin)( 
    int Mb, int Nb, int k, int l )
{
  TSp_mat<DOUBLE> *A = new TSp_mat<DOUBLE>(Mb*k, Nb*l);

  TSp_MAT_SET_DOUBLE(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_DOUBLE_NAME(uscr_variable_block_begin)( 
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<DOUBLE> *A = new TSp_mat<DOUBLE>(  
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_DOUBLE(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  DOUBLE Level 2/3 insertion routines */

int BLAS_DOUBLE_NAME(uscr_insert_entry)( 
    blas_sparse_matrix A, SPBLAS_DOUBLE_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_DOUBLE_IN( val ), i, j);
}

int BLAS_DOUBLE_NAME(uscr_insert_entries)( 
    blas_sparse_matrix A, int nz, 
    SPBLAS_VECTOR_DOUBLE_IN val, 
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_DOUBLE_IN( val ), indx, jndx);
}

int BLAS_DOUBLE_NAME(uscr_insert_col)( 
    blas_sparse_matrix A, int j, int nz, 
    SPBLAS_VECTOR_DOUBLE_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_DOUBLE_IN( val ), indx);
}

int BLAS_DOUBLE_NAME(uscr_insert_row)( 
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_DOUBLE_IN val, const int *indx );

int BLAS_DOUBLE_NAME(uscr_insert_clique)( 
    blas_sparse_matrix A, 
    const int k, 
    const int l,
    SPBLAS_VECTOR_DOUBLE_IN val, 
    const int row_stride,
    const int col_stride, 
    const int *indx,
    const int *jndx );

int BLAS_DOUBLE_NAME(uscr_insert_block)( 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_DOUBLE_IN val,
    int row_stride, 
    int col_stride, 
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_DOUBLE_IN( val ), 
        row_stride, col_stride, i, j);
}

int BLAS_DOUBLE_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  DOUBLE Level 2/3 computational routines */

 int BLAS_DOUBLE_NAME(usmv)(enum 
    blas_trans_type transa, 
    SPBLAS_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_DOUBLE_IN x, 
    int incx, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT y, 
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_DOUBLE_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_DOUBLE_IN( x ), incx, 
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( y ), incy);
}

int BLAS_DOUBLE_NAME(usmm)( 
    enum blas_order_type order, 
    enum blas_trans_type transa,
    int nrhs, 
    SPBLAS_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_DOUBLE_IN b,
    int ldb, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT c, 
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs, 
      SPBLAS_TO_DOUBLE_IN( alpha), A, 
      SPBLAS_TO_VECTOR_DOUBLE_IN(b), ldb, 
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( c ), ldc);
}

int BLAS_DOUBLE_NAME(ussv)( 
    enum blas_trans_type transa, 
    SPBLAS_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT x, 
    int incx )
{
  return BLAS_xussv( transa, 
        SPBLAS_TO_DOUBLE_IN( alpha ), A, 
        SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( x ), 
        incx);
}


int BLAS_DOUBLE_NAME(ussm)( 
    enum blas_order_type order, 
    enum blas_trans_type transt,
    int nrhs, 
    SPBLAS_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_DOUBLE_IN_OUT b, 
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs, 
      SPBLAS_TO_DOUBLE_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_DOUBLE_IN_OUT( b ), ldb);
}


/*  ----   end of DOUBLE routines -------  */




 void BLAS_COMPLEX_DOUBLE_NAME(usdot)( 
    enum blas_conj_type conj_flag, 
    int  nz, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x,  
    const int *index,  
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN y, 
    int incy, 
    SPBLAS_COMPLEX_DOUBLE_OUT r, 
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz, 
          SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ), index, 
          SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( y ), incy, 
          SPBLAS_TO_COMPLEX_DOUBLE_OUT( r ), index_base);
}

 void BLAS_COMPLEX_DOUBLE_NAME(usaxpy)(
    int nz, 
    SPBLAS_COMPLEX_DOUBLE_IN alpha, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x, 
    const int *index, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y, 
    int incy, 
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ), 
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ), index, 
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( y ), 
      incy, index_base);
}

 void BLAS_COMPLEX_DOUBLE_NAME(usga)( 
    int nz, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN y, 
    int incy, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( y ), incy, 
        SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( x ), indx, index_base);

}

 void BLAS_COMPLEX_DOUBLE_NAME(usgz)( 
    int nz, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y, 
    int incy, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_COMPLEX_DOUBLE_IN_OUT(y ), incy, SPBLAS_TO_COMPLEX_DOUBLE_OUT(  x ), 
    indx, index_base);
}


 void BLAS_COMPLEX_DOUBLE_NAME(ussc)(
    int nz, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y, 
    int incy, 
    const int *index, 
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ), 
  SPBLAS_TO_COMPLEX_DOUBLE_IN_OUT( y ), incy, index, 
  index_base);
}

/*  COMPLEX_DOUBLE Level 2/3 creation routines */

int BLAS_COMPLEX_DOUBLE_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<COMPLEX_DOUBLE> *A = new TSp_mat<COMPLEX_DOUBLE>(M, N);
  TSp_MAT_SET_COMPLEX_DOUBLE(A);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_DOUBLE_NAME(uscr_block_begin)( 
    int Mb, int Nb, int k, int l )
{
  TSp_mat<COMPLEX_DOUBLE> *A = new TSp_mat<COMPLEX_DOUBLE>(Mb*k, Nb*l);

  TSp_MAT_SET_COMPLEX_DOUBLE(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_DOUBLE_NAME(uscr_variable_block_begin)( 
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<COMPLEX_DOUBLE> *A = new TSp_mat<COMPLEX_DOUBLE>(  
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_COMPLEX_DOUBLE(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  COMPLEX_DOUBLE Level 2/3 insertion routines */

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_entry)( 
    blas_sparse_matrix A, SPBLAS_COMPLEX_DOUBLE_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_COMPLEX_DOUBLE_IN( val ), i, j);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_entries)( 
    blas_sparse_matrix A, int nz, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val, 
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( val ), indx, jndx);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_col)( 
    blas_sparse_matrix A, int j, int nz, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( val ), indx);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_row)( 
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val, const int *indx );

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_clique)( 
    blas_sparse_matrix A, 
    const int k, 
    const int l,
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val, 
    const int row_stride,
    const int col_stride, 
    const int *indx,
    const int *jndx );

int BLAS_COMPLEX_DOUBLE_NAME(uscr_insert_block)( 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN val,
    int row_stride, 
    int col_stride, 
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( val ), 
        row_stride, col_stride, i, j);
}

int BLAS_COMPLEX_DOUBLE_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  COMPLEX_DOUBLE Level 2/3 computational routines */

 int BLAS_COMPLEX_DOUBLE_NAME(usmv)(enum 
    blas_trans_type transa, 
    SPBLAS_COMPLEX_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN x, 
    int incx, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT y, 
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN( x ), incx, 
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( y ), incy);
}

int BLAS_COMPLEX_DOUBLE_NAME(usmm)( 
    enum blas_order_type order, 
    enum blas_trans_type transa,
    int nrhs, 
    SPBLAS_COMPLEX_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN b,
    int ldb, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT c, 
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs, 
      SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha), A, 
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN(b), ldb, 
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( c ), ldc);
}

int BLAS_COMPLEX_DOUBLE_NAME(ussv)( 
    enum blas_trans_type transa, 
    SPBLAS_COMPLEX_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT x, 
    int incx )
{
  return BLAS_xussv( transa, 
        SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ), A, 
        SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( x ), 
        incx);
}


int BLAS_COMPLEX_DOUBLE_NAME(ussm)( 
    enum blas_order_type order, 
    enum blas_trans_type transt,
    int nrhs, 
    SPBLAS_COMPLEX_DOUBLE_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_DOUBLE_IN_OUT b, 
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs, 
      SPBLAS_TO_COMPLEX_DOUBLE_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_COMPLEX_DOUBLE_IN_OUT( b ), ldb);
}




/*  ----   end of COMPLEX_COMPLEX_COMPLEX_DOUBLE routines -------  */


/*********/
/*  ----  double  Level 1 rouintes ----- */
/*********/

 void BLAS_FLOAT_NAME(usdot)( 
    enum blas_conj_type conj_flag, 
    int  nz, 
    SPBLAS_VECTOR_FLOAT_IN x,  
    const int *index,  
    SPBLAS_VECTOR_FLOAT_IN y, 
    int incy, 
    SPBLAS_FLOAT_OUT r, 
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz, 
          SPBLAS_TO_VECTOR_FLOAT_IN( x ), index, 
          SPBLAS_TO_VECTOR_FLOAT_IN( y ), incy, 
          SPBLAS_TO_FLOAT_OUT( r ), index_base);
}

 void BLAS_FLOAT_NAME(usaxpy)(
    int nz, 
    SPBLAS_FLOAT_IN alpha, 
    SPBLAS_VECTOR_FLOAT_IN x, 
    const int *index, 
    SPBLAS_VECTOR_FLOAT_IN_OUT y, 
    int incy, 
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_FLOAT_IN( alpha ), 
      SPBLAS_TO_VECTOR_FLOAT_IN( x ), index, 
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( y ), 
      incy, index_base);
}

 void BLAS_FLOAT_NAME(usga)( 
    int nz, 
    SPBLAS_VECTOR_FLOAT_IN y, 
    int incy, 
    SPBLAS_VECTOR_FLOAT_IN_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_FLOAT_IN( y ), incy, 
        SPBLAS_TO_VECTOR_FLOAT_IN_OUT( x ), indx, index_base);

}

 void BLAS_FLOAT_NAME(usgz)( 
    int nz, 
    SPBLAS_VECTOR_FLOAT_IN_OUT y, 
    int incy, 
    SPBLAS_VECTOR_FLOAT_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_FLOAT_IN_OUT(y ), incy, SPBLAS_TO_FLOAT_OUT(  x ), 
    indx, index_base);
}


 void BLAS_FLOAT_NAME(ussc)(
    int nz, 
    SPBLAS_VECTOR_FLOAT_IN x, 
    SPBLAS_VECTOR_FLOAT_IN_OUT y, 
    int incy, 
    const int *index, 
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_FLOAT_IN( x ), 
  SPBLAS_TO_FLOAT_IN_OUT( y ), incy, index, 
  index_base);
}

/*  FLOAT Level 2/3 creation routines */

int BLAS_FLOAT_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<FLOAT> *A = new TSp_mat<FLOAT>(M, N);
  TSp_MAT_SET_FLOAT(A);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_FLOAT_NAME(uscr_block_begin)( 
    int Mb, int Nb, int k, int l )
{
  TSp_mat<FLOAT> *A = new TSp_mat<FLOAT>(Mb*k, Nb*l);

  TSp_MAT_SET_FLOAT(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_FLOAT_NAME(uscr_variable_block_begin)( 
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<FLOAT> *A = new TSp_mat<FLOAT>(  
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_FLOAT(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  FLOAT Level 2/3 insertion routines */

int BLAS_FLOAT_NAME(uscr_insert_entry)( 
    blas_sparse_matrix A, SPBLAS_FLOAT_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_FLOAT_IN( val ), i, j);
}

int BLAS_FLOAT_NAME(uscr_insert_entries)( 
    blas_sparse_matrix A, int nz, 
    SPBLAS_VECTOR_FLOAT_IN val, 
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_FLOAT_IN( val ), indx, jndx);
}

int BLAS_FLOAT_NAME(uscr_insert_col)( 
    blas_sparse_matrix A, int j, int nz, 
    SPBLAS_VECTOR_FLOAT_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_FLOAT_IN( val ), indx);
}

int BLAS_FLOAT_NAME(uscr_insert_row)( 
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_FLOAT_IN val, const int *indx );

int BLAS_FLOAT_NAME(uscr_insert_clique)( 
    blas_sparse_matrix A, 
    const int k, 
    const int l,
    SPBLAS_VECTOR_FLOAT_IN val, 
    const int row_stride,
    const int col_stride, 
    const int *indx,
    const int *jndx );

int BLAS_FLOAT_NAME(uscr_insert_block)( 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_FLOAT_IN val,
    int row_stride, 
    int col_stride, 
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_FLOAT_IN( val ), 
        row_stride, col_stride, i, j);
}

int BLAS_FLOAT_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  FLOAT Level 2/3 computational routines */

 int BLAS_FLOAT_NAME(usmv)(enum 
    blas_trans_type transa, 
    SPBLAS_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_FLOAT_IN x, 
    int incx, 
    SPBLAS_VECTOR_FLOAT_IN_OUT y, 
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_FLOAT_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_FLOAT_IN( x ), incx, 
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( y ), incy);
}

int BLAS_FLOAT_NAME(usmm)( 
    enum blas_order_type order, 
    enum blas_trans_type transa,
    int nrhs, 
    SPBLAS_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_FLOAT_IN b,
    int ldb, 
    SPBLAS_VECTOR_FLOAT_IN_OUT c, 
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs, 
      SPBLAS_TO_FLOAT_IN( alpha), A, 
      SPBLAS_TO_VECTOR_FLOAT_IN(b), ldb, 
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( c ), ldc);
}

int BLAS_FLOAT_NAME(ussv)( 
    enum blas_trans_type transa, 
    SPBLAS_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_FLOAT_IN_OUT x, 
    int incx )
{
  return BLAS_xussv( transa, 
        SPBLAS_TO_FLOAT_IN( alpha ), A, 
        SPBLAS_TO_VECTOR_FLOAT_IN_OUT( x ), 
        incx);
}


int BLAS_FLOAT_NAME(ussm)( 
    enum blas_order_type order, 
    enum blas_trans_type transt,
    int nrhs, 
    SPBLAS_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_FLOAT_IN_OUT b, 
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs, 
      SPBLAS_TO_FLOAT_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_FLOAT_IN_OUT( b ), ldb);
}


/*  ----   end of FLOAT routines -------  */




 void BLAS_COMPLEX_FLOAT_NAME(usdot)( 
    enum blas_conj_type conj_flag, 
    int  nz, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x,  
    const int *index,  
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN y, 
    int incy, 
    SPBLAS_COMPLEX_FLOAT_OUT r, 
    enum blas_base_type index_base)
{
   BLAS_xusdot(conj_flag, nz, 
          SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ), index, 
          SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( y ), incy, 
          SPBLAS_TO_COMPLEX_FLOAT_OUT( r ), index_base);
}

 void BLAS_COMPLEX_FLOAT_NAME(usaxpy)(
    int nz, 
    SPBLAS_COMPLEX_FLOAT_IN alpha, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x, 
    const int *index, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y, 
    int incy, 
    enum blas_base_type index_base)
{
  BLAS_xusaxpy(nz, SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ), 
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ), index, 
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( y ), 
      incy, index_base);
}

 void BLAS_COMPLEX_FLOAT_NAME(usga)( 
    int nz, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN y, 
    int incy, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusga( nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( y ), incy, 
        SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( x ), indx, index_base);

}

 void BLAS_COMPLEX_FLOAT_NAME(usgz)( 
    int nz, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y, 
    int incy, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_OUT x, 
    const int *indx,
    enum blas_base_type index_base )
{
  BLAS_xusgz(nz, SPBLAS_TO_COMPLEX_FLOAT_IN_OUT(y ), incy, SPBLAS_TO_COMPLEX_FLOAT_OUT(  x ), 
    indx, index_base);
}


 void BLAS_COMPLEX_FLOAT_NAME(ussc)(
    int nz, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y, 
    int incy, 
    const int *index, 
    enum blas_base_type index_base)
{
  BLAS_xussc(nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ), 
  SPBLAS_TO_COMPLEX_FLOAT_IN_OUT( y ), incy, index, 
  index_base);
}

/*  COMPLEX_FLOAT Level 2/3 creation routines */

int BLAS_COMPLEX_FLOAT_NAME(uscr_begin)(int M, int N)
{
  TSp_mat<COMPLEX_FLOAT> *A = new TSp_mat<COMPLEX_FLOAT>(M, N);
  TSp_MAT_SET_COMPLEX_FLOAT(A);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_FLOAT_NAME(uscr_block_begin)( 
    int Mb, int Nb, int k, int l )
{
  TSp_mat<COMPLEX_FLOAT> *A = new TSp_mat<COMPLEX_FLOAT>(Mb*k, Nb*l);

  TSp_MAT_SET_COMPLEX_FLOAT(A);
  A->set_const_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);
}

blas_sparse_matrix BLAS_COMPLEX_FLOAT_NAME(uscr_variable_block_begin)( 
    int Mb, int Nb, const int *k, const int *l )
{
  TSp_mat<COMPLEX_FLOAT> *A = new TSp_mat<COMPLEX_FLOAT>(  
                std::accumulate(k, k+Mb, 0), std::accumulate(l, l+Nb, 0) );

  TSp_MAT_SET_COMPLEX_FLOAT(A);
  A->set_var_block_parameters(Mb, Nb, k, l);

  return Table_insert(A);

}

/*  COMPLEX_FLOAT Level 2/3 insertion routines */

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_entry)( 
    blas_sparse_matrix A, SPBLAS_COMPLEX_FLOAT_IN val, int i, int j )
{
  return BLAS_xuscr_insert_entry(A, SPBLAS_TO_COMPLEX_FLOAT_IN( val ), i, j);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_entries)( 
    blas_sparse_matrix A, int nz, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val, 
    const int *indx, const int *jndx )
{
  return BLAS_xuscr_insert_entries(A, nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( val ), indx, jndx);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_col)( 
    blas_sparse_matrix A, int j, int nz, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val, const int *indx )
{
  return BLAS_xuscr_insert_col(A, j, nz, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( val ), indx);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_row)( 
  blas_sparse_matrix A, int i, int nz,
  SPBLAS_VECTOR_COMPLEX_FLOAT_IN val, const int *indx );

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_clique)( 
    blas_sparse_matrix A, 
    const int k, 
    const int l,
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val, 
    const int row_stride,
    const int col_stride, 
    const int *indx,
    const int *jndx );

int BLAS_COMPLEX_FLOAT_NAME(uscr_insert_block)( 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN val,
    int row_stride, 
    int col_stride, 
    int i, int j )
{
  return BLAS_xuscr_insert_block(
        A, SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( val ), 
        row_stride, col_stride, i, j);
}

int BLAS_COMPLEX_FLOAT_NAME(uscr_end)(blas_sparse_matrix A)
{
  return BLAS_xuscr_end(A);
}

/*  COMPLEX_FLOAT Level 2/3 computational routines */

 int BLAS_COMPLEX_FLOAT_NAME(usmv)(enum 
    blas_trans_type transa, 
    SPBLAS_COMPLEX_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN x, 
    int incx, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT y, 
    int incy )
{
  return BLAS_xusmv(
      transa,   SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN( x ), incx, 
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( y ), incy);
}

int BLAS_COMPLEX_FLOAT_NAME(usmm)( 
    enum blas_order_type order, 
    enum blas_trans_type transa,
    int nrhs, 
    SPBLAS_COMPLEX_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN b,
    int ldb, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT c, 
    int ldc )
{
  return BLAS_xusmm(
      order, transa, nrhs, 
      SPBLAS_TO_COMPLEX_FLOAT_IN( alpha), A, 
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN(b), ldb, 
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( c ), ldc);
}

int BLAS_COMPLEX_FLOAT_NAME(ussv)( 
    enum blas_trans_type transa, 
    SPBLAS_COMPLEX_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT x, 
    int incx )
{
  return BLAS_xussv( transa, 
        SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ), A, 
        SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( x ), 
        incx);
}


int BLAS_COMPLEX_FLOAT_NAME(ussm)( 
    enum blas_order_type order, 
    enum blas_trans_type transt,
    int nrhs, 
    SPBLAS_COMPLEX_FLOAT_IN alpha, 
    blas_sparse_matrix A, 
    SPBLAS_VECTOR_COMPLEX_FLOAT_IN_OUT b, 
    int ldb )
{
  return BLAS_xussm(order, transt, nrhs, 
      SPBLAS_TO_COMPLEX_FLOAT_IN( alpha ), A, 
      SPBLAS_TO_VECTOR_COMPLEX_FLOAT_IN_OUT( b ), ldb);
}




/*  ----   end of COMPLEX_COMPLEX_COMPLEX_FLOAT routines -------  */



















