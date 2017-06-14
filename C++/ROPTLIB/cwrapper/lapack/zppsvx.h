#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zppsvx_(char *fact, char *uplo, integer *n, integer *nrhs, doublecomplex *ap, doublecomplex *afp, char *equed, doublereal *s, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#ifdef __cplusplus
}
#endif