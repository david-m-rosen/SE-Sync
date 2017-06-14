#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n, complex *a, integer *lda, complex *b, integer *ldb, integer *sdim, complex *alpha, complex *beta, complex *vsl, integer *ldvsl, complex *vsr, integer *ldvsr, real *rconde, real *rcondv, complex *work, integer *lwork, real *rwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif