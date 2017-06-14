#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int sgeesx_(char *jobvs, char *sort, L_fp select, char *sense, integer *n, real *a, integer *lda, integer *sdim, real *wr, real *wi, real *vs, integer *ldvs, real *rconde, real *rcondv, real *work, integer *lwork, integer *iwork, integer *liwork, logical *bwork, integer *info);

#ifdef __cplusplus
}
#endif