#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cggev_(char *jobvl, char *jobvr, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *vl, integer *ldvl, complex *vr, integer *ldvr, complex *work, integer *lwork, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif