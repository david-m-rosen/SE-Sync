#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int ctbrfs_(char *uplo, char *trans, char *diag, integer *n, integer *kd, integer *nrhs, complex *ab, integer *ldab, complex *b, integer *ldb, complex *x, integer *ldx, real *ferr, real *berr, complex *work, real *rwork, integer *info);

#ifdef __cplusplus
}
#endif