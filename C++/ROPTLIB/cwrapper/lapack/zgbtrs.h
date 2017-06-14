#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgbtrs_(char *trans, integer *n, integer *kl, integer *ku, integer *nrhs, doublecomplex *ab, integer *ldab, integer *ipiv, doublecomplex *b, integer *ldb, integer *info);

#ifdef __cplusplus
}
#endif