#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cggevx_(char *balanc, char *jobvl, char *jobvr, char *sense, integer *n, complex *a, integer *lda, complex *b, integer *ldb, complex *alpha, complex *beta, complex *vl, integer *ldvl, complex *vr, integer *ldvr, integer *ilo, integer *ihi, real *lscale, real *rscale, real *abnrm, real *bbnrm, real *rconde, real *rcondv, complex *work, integer *lwork, real *rwork, integer *iwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif