#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cgges_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, integer *n, complex *a, integer *lda, complex *b, integer *ldb, integer *sdim, complex *alpha, complex *beta, complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, complex *work, integer *lwork, real *rwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif