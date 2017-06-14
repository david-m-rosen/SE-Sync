#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgesvx_(char *fact, char *trans, integer *n, integer *nrhs, real *a, integer *lda, real *af, integer *ldaf, integer *ipiv, char *equed, real *r__, real *c__, real *b, integer *ldb, real *x, integer *ldx, real *rcond, real *ferr, real *berr, real *work, integer *iwork, integer *info);

#ifdef __cplusplus
}
#endif