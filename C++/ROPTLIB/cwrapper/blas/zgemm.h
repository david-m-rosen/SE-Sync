#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int zgemm_(char *transa, char *transb, integer *m, integer *n, integer *k, doublecomplex *alpha, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, doublecomplex *beta, doublecomplex *c__, integer *ldc);

#ifdef __cplusplus
}
#endif