#ifndef BLAS_SPARSE_PROTO_H
#define BLAS_SPARSE_PROTO_H

typedef int blas_sparse_matrix;


  /* Level 1 Computational Routines */

void BLAS_susdot( enum blas_conj_type conj, int nz, const float *x, 
                  const int *indx, const float *y, int incy, float *r,
                  enum blas_base_type index_base );
void BLAS_dusdot( enum blas_conj_type conj, int nz, const double *x, 
                  const int *indx, const double *y, int incy, double *r,
                  enum blas_base_type index_base );
void BLAS_cusdot( enum blas_conj_type conj, int nz, const void *x, 
                  const int *indx, const void *y, int incy, void *r,
                  enum blas_base_type index_base );
void BLAS_zusdot( enum blas_conj_type conj, int nz, const void *x, 
                  const int *indx, const void *y, int incy, void *r,
                  enum blas_base_type index_base );

void BLAS_susaxpy( int nz, float alpha, const float *x, const int *indx,
                 float *y, int incy, enum blas_base_type index_base );
void BLAS_dusaxpy( int nz, double alpha, const double *x, const int *indx,
                 double *y, int incy, enum blas_base_type index_base );
void BLAS_cusaxpy( int nz, const void *alpha, const void *x, const int *indx,
                 void *y, int incy, enum blas_base_type index_base );
void BLAS_zusaxpy( int nz, const void *alpha, const void *x, const int *indx,
                 void *y, int incy, enum blas_base_type index_base );

void BLAS_susga( int nz, const float *y, int incy, float *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_dusga( int nz, const double *y, int incy, double *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_cusga( int nz, const void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_zusga( int nz, const void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );

void BLAS_susgz( int nz, float *y, int incy, float *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_dusgz( int nz, double *y, int incy, double *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_cusgz( int nz, void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );
void BLAS_zusgz( int nz, void *y, int incy, void *x, const int *indx,
              enum blas_base_type index_base );

void BLAS_sussc( int nz, const float *x, float *y, int incy, const int *indx,
              enum blas_base_type index_base );
void BLAS_dussc( int nz, const double *x, double *y, int incy, const int *indx,
              enum blas_base_type index_base );
void BLAS_cussc( int nz, const void *x, void *y, int incy, const int *indx,
              enum blas_base_type index_base );
void BLAS_zussc( int nz, const void *x, void *y, int incy, const int *indx,
              enum blas_base_type index_base );

               /* Level 2 Computational Routines */

int BLAS_susmv( enum blas_trans_type transa, float alpha, 
    blas_sparse_matrix A, const float *x, int incx, float *y, int incy );
int BLAS_dusmv( enum blas_trans_type transa, double alpha, 
    blas_sparse_matrix A, const double *x, int incx, double *y, int incy );
int BLAS_cusmv( enum blas_trans_type transa, const void *alpha, 
    blas_sparse_matrix A, const void *x, int incx, void *y, int incy );
int BLAS_zusmv( enum blas_trans_type transa, const void *alpha, 
    blas_sparse_matrix A, const void *x, int incx, void *y, int incy );

int BLAS_sussv( enum blas_trans_type transt, float alpha, 
    blas_sparse_matrix T, float *x, int incx );
int BLAS_dussv( enum blas_trans_type transt, double alpha, 
    blas_sparse_matrix T, double *x, int incx );
int BLAS_cussv( enum blas_trans_type transt, const void *alpha, 
    blas_sparse_matrix T, void *x, int incx );
int BLAS_zussv( enum blas_trans_type transt, const void *alpha, 
    blas_sparse_matrix T, void *x, int incx );

               /* Level 3 Computational Routines */

int BLAS_susmm( enum blas_order_type order, enum blas_trans_type transa,
    int nrhs, float alpha, blas_sparse_matrix A, const float *b, int ldb,
        float *c, int ldc );
int BLAS_dusmm( enum blas_order_type order, enum blas_trans_type transa,
        int nrhs, double alpha, blas_sparse_matrix A, const double *b,
        int ldb, double *c, int ldc );
int BLAS_cusmm( enum blas_order_type order, enum blas_trans_type transa,
         int nrhs, const void *alpha, blas_sparse_matrix A, const void *b, 
     int ldb, void *c, int ldc );
int BLAS_zusmm( enum blas_order_type order, enum blas_trans_type transa,
         int nrhs, const void *alpha, blas_sparse_matrix A, const void *b, 
     int ldb, void *c, int ldc );

int BLAS_sussm( enum blas_order_type order, enum blas_trans_type transt,
              int nrhs, float alpha, int t, float *b, int ldb );
int BLAS_dussm( enum blas_order_type order, enum blas_trans_type transt,
              int nrhs, double alpha, int t, double *b, int ldb );
int BLAS_cussm( enum blas_order_type order, enum blas_trans_type transt,
              int nrhs, const void *alpha, int t, void *b, int ldb );
int BLAS_zussm( enum blas_order_type order, enum blas_trans_type transt,
              int nrhs, const void *alpha, int t, void *b, int ldb );

               /* Handle Management Routines */

               /* Creation Routines */

blas_sparse_matrix BLAS_suscr_begin( int m, int n );
blas_sparse_matrix BLAS_duscr_begin( int m, int n );
blas_sparse_matrix BLAS_cuscr_begin( int m, int n );
blas_sparse_matrix BLAS_zuscr_begin( int m, int n );


blas_sparse_matrix BLAS_suscr_block_begin( int Mb, int Nb, int k, int l );
blas_sparse_matrix BLAS_duscr_block_begin( int Mb, int Nb, int k, int l );
blas_sparse_matrix BLAS_cuscr_block_begin( int Mb, int Nb, int k, int l );
blas_sparse_matrix BLAS_zuscr_block_begin( int Mb, int Nb, int k, int l );

blas_sparse_matrix BLAS_suscr_variable_block_begin( int Mb, int Nb, 
		const int *k, const int *l );
blas_sparse_matrix BLAS_duscr_variable_block_begin( int Mb, int Nb, 
		const int *k, const int *l );
blas_sparse_matrix BLAS_cuscr_variable_block_begin( int Mb, int Nb, 
		const int *k, const int *l );
blas_sparse_matrix BLAS_zuscr_variable_block_begin( int Mb, int Nb, 
		const int *k, const int *l );


               /* Insertion Routines */

int BLAS_suscr_insert_entry( blas_sparse_matrix A, float val, int i, int j );
int BLAS_duscr_insert_entry( blas_sparse_matrix A, double val, int i, int j );
int BLAS_cuscr_insert_entry( blas_sparse_matrix A, const void *val, int i, int j );
int BLAS_zuscr_insert_entry( blas_sparse_matrix A, const void *val, int i, int j );

int BLAS_suscr_insert_entries( blas_sparse_matrix A, int nz, const float *val,
                            const int *indx, const int *jndx );
int BLAS_duscr_insert_entries( blas_sparse_matrix A, int nz, const double *val,
                            const int *indx, const int *jndx );
int BLAS_cuscr_insert_entries( blas_sparse_matrix A, int nz, const void *val,
                            const int *indx, const int *jndx );
int BLAS_zuscr_insert_entries( blas_sparse_matrix A, int nz, const void *val,
                            const int *indx, const int *jndx );

int BLAS_suscr_insert_col( blas_sparse_matrix A, int j, int nz,
                           const float *val, const int *indx );
int BLAS_duscr_insert_col( blas_sparse_matrix A, int j, int nz,
                           const double *val, const int *indx );
int BLAS_cuscr_insert_col( blas_sparse_matrix A, int j, int nz,
                           const void *val, const int *indx );
int BLAS_zuscr_insert_col( blas_sparse_matrix A, int j, int nz,
                           const void *val, const int *indx );

int BLAS_suscr_insert_row( blas_sparse_matrix A, int i, int nz,
                           const float *val, const int *indx );
int BLAS_duscr_insert_row( blas_sparse_matrix A, int i, int nz,
                           const double *val, const int *indx );
int BLAS_cuscr_insert_row( blas_sparse_matrix A, int i, int nz,
                           const void *val, const int *indx );
int BLAS_zuscr_insert_row( blas_sparse_matrix A, int i, int nz,
                           const void *val, const int *indx );

int BLAS_suscr_insert_clique( blas_sparse_matrix A, const int k, const int l, 
                        const float *val, const int row_stride, 
                        const int col_stride, const int *indx, 
                        const int *jndx );
int BLAS_duscr_insert_clique( blas_sparse_matrix A, const int k, const int l, 
                        const double *val, const int row_stride, 
                        const int col_stride, const int *indx, 
                        const int *jndx );
int BLAS_cuscr_insert_clique( blas_sparse_matrix A, const int k, const int l, 
                        const void *val, const int row_stride, 
                        const int col_stride, const int *indx, 
                        const int *jndx );
int BLAS_zuscr_insert_clique( blas_sparse_matrix A, const int k, const int l, 
                        const void *val, const int row_stride, 
                        const int col_stride, const int *indx, 
                        const int *jndx );

int BLAS_suscr_insert_block( blas_sparse_matrix A, const float *val, 
                        int row_stride, int col_stride, int i, int j );
int BLAS_duscr_insert_block( blas_sparse_matrix A, const double *val, 
                        int row_stride, int col_stride, int i, int j );
int BLAS_cuscr_insert_block( blas_sparse_matrix A, const void *val, 
                        int row_stride, int col_stride, int i, int j );
int BLAS_zuscr_insert_block( blas_sparse_matrix A, const void *val, 
                        int row_stride, int col_stride, int i, int j );

               /* Completion of Construction Routines */

int BLAS_suscr_end( blas_sparse_matrix A );
int BLAS_duscr_end( blas_sparse_matrix A );
int BLAS_cuscr_end( blas_sparse_matrix A );
int BLAS_zuscr_end( blas_sparse_matrix A );

               /* Matrix Property Routines */

int BLAS_usgp( blas_sparse_matrix A, int pname );

int BLAS_ussp( blas_sparse_matrix A, int pname );

               /* Destruction Routine */

int BLAS_usds( blas_sparse_matrix A );

#endif
  /* BLAS_SPARSE_PROTO_H */
