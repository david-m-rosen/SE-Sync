#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sggesx_(char *jobvsl, char *jobvsr, char *sort, L_fp selctg, char *sense, integer *n, real *a, integer *lda, real *b, integer *ldb, integer *sdim, real *alphar, real *alphai, real *beta, real *vsl, integer *ldvsl, real *vsr, integer *ldvsr, real *rconde, real *rcondv, real *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif