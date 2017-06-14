#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zpbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, doublecomplex *ab, integer *ldab, doublecomplex *afb, integer *ldafb, char *equed, doublereal *s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#ifdef __cplusplus
}
#endif