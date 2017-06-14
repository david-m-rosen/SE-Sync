#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zpprfs_(char *uplo, integer *n, integer *nrhs, doublecomplex *ap, doublecomplex *afp, doublecomplex *b, integer *ldb, doublecomplex *x, integer *ldx, doublereal *ferr, doublereal *berr, doublecomplex *work, doublereal *rwork, integer *info);

#ifdef __cplusplus
}
#endif