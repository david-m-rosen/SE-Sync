#ifdef __cplusplus
extern "C" { 
#endif  

#include "f2c.h" 

int chbgvd_(char *jobz, char *uplo, integer *n, integer *ka, integer *kb, complex *ab, integer *ldab, complex *bb, integer *ldbb, real *w, complex *z__, integer *ldz, complex *work, integer *lwork, real *rwork, integer *lrwork, integer *iwork, integer *liwork, integer *info);

#ifdef __cplusplus
}
#endif