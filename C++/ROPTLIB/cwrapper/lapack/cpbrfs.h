#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int cpbrfs_(char *uplo, integer *n, integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *afb, integer *ldafb, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif