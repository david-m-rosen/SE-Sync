#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgees_(char *jobvs, char *sort, L_fp select, integer *n, doublecomplex *a, integer *lda, integer *sdim, doublecomplex *w, doublecomplex *vs, integer *ldvs, doublecomplex *work, integer *lwork, doublereal *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif