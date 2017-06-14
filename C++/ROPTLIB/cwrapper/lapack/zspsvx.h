#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zspsvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublecomplex *ap, doublecomplex *afp, integer *ipiv, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#ifdef __cplusplus
}
#endif