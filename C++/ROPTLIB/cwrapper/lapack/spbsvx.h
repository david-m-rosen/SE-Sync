#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int spbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, real *ab, integer *ldab, real *afb, integer *ldafb, char *equed, real *s, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif