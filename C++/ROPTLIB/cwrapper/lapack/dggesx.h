#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int dggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp delctg, char *sense, integer *n, doublereal *a, integer *lda, doublereal *b, integer *ldb, integer *sdim, doublereal *alphar, doublereal *alphai, doublereal *beta, doublereal *vsl, integer *ldvsl, doublereal *vsr, integer *ldvsr, doublereal *rconde, doublereal *rcondv, doublereal *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif