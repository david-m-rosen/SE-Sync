#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dpbsvx_(char *fact, char *uplo, integer *n, integer *kd, integer *nrhs, doublereal *ab, integer *ldab, doublereal *afb, integer *ldafb, char *equed, doublereal *s, doublereal *b, integer *ldb, doublereal *x, integer *ldx, doublereal *rcond, doublereal *ferr, doublereal *berr, doublereal *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif