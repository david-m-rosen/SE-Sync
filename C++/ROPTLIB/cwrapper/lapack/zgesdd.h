#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgesdd_(char *jobz, integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, integer *lwork, doublereal *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif