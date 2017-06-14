#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp delctg, char *sense, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *sdim, doublecomplex *alpha, doublecomplex *beta, doublecomplex *vsl, integer *ldvsl, doublecomplex *vsr, integer *ldvsr, doublereal *rconde, doublereal *rcondv, doublecomplex *work, integer *lwork, doublereal *rwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif