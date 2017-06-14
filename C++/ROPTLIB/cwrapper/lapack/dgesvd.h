#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dgesvd_(char *jobu, char *jobvt, integer *m, integer *n, doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, integer *info);

#ifdef __cplusplus
}
#endif