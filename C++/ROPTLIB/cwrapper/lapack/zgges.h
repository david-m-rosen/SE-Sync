#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgges_(char *jobvsl, char *jobvsr, char *sort, L_fp delctg, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif