#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublecomplex *a, integer *lda, doublereal *s, doublecomplex *u, integer *ldu, doublecomplex *vt, integer *ldvt, doublecomplex *work, integer *lwork, doublereal *rwork, integer *info);

#ifdef __cplusplus
}
#endif