#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zhegvx_(integer *itype, char *jobz, char *range, char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublereal *vl, doublereal *vu, integer *il, integer *iu, doublereal *abstol, integer *m, doublereal *w, doublecomplex *z__, integer *ldz, doublecomplex *work, integer *lwork, doublereal *rwork, integer *iwork, integer *ifail, integer *info);

#ifdef __cplusplus
}
#endif