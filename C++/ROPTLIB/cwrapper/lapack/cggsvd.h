#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cggsvd_(char *jobu, char *jobv, char *jobq, integer *m, integer *n, integer *p, integer *k, integer *l, complex *a, integer *lda, complex *b, integer *ldb, real *alpha, real *beta, complex *u, integer *ldu, complex *v, integer *ldv, complex *q, integer *ldq, complex *work, real *rwork, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif