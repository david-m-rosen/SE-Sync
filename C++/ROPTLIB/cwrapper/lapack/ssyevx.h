#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ssyevx_(char *jobz, char *range, char *uplo, integer *n, real *a, integer *lda, real *vl, real *vu, integer *il, integer *iu, real *abstol, integer *m, real *w, real *z__, integer *ldz, real *work, integer *lwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif