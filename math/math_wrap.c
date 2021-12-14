// external C prototypes for BLAS & LAPACK routines used by MOPAC
#include <stddef.h>

void daxpy_(int*, double*, double*, int*, double*, int*);
void dcopy_(int*, double*, int*, double*, int*);
double ddot_(int*, double*, int*, double*, int*);
void dgemm_(char*, char*, int*, int*, int*, double*, double*, int*, double* ,int* ,double*, double* ,int*);
void dgemv_(char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
void dgesv_(int*, int*, double*, int*, int*, double*, int*, int*);
void dlacpy_(char*, int*, int*, double*, int*, double*, int*);
void dpotrf_(char*, int*, double*, int*, int*);
void dpotri_(char*, int*, double*, int*, int*);
void drot_(int*, double*, int*, double*, int*, double*, double*);
void dscal_(int*, double*, double*, int*);
void dswap_(int*, double*, int*, double*, int*);
void dsyevd_(char*, char*, int*, double*, int*, double*, double*, int*, int*, int* ,int*);
void dsymm_(char*, char*, int*, int*, double*, double*, int*, double*, int*, double*, double*, int*);
void dsyrk_(char*, char*, int*, int*, double*, double*, int*, double*, double*, int*);
void dtpttr_(char*, int*, double*, double*, int*, int*);
void dtrttp_(char*, int*, double*, int*, double*, int*);
int idamax_(int*, double*, int*);

// touch the functions to force their linkage
void touch()
{
  daxpy_(NULL, NULL, NULL, NULL, NULL, NULL);
  dcopy_(NULL, NULL, NULL, NULL, NULL);
  ddot_(NULL, NULL, NULL, NULL, NULL);
  dgemm_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dgemv_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dgesv_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dlacpy_(NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dpotrf_(NULL, NULL, NULL, NULL, NULL);
  dpotri_(NULL, NULL, NULL, NULL, NULL);
  drot_(NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dscal_(NULL, NULL, NULL, NULL);
  dswap_(NULL, NULL, NULL, NULL, NULL);
  dsyevd_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dsymm_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dsyrk_(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL);
  dtpttr_(NULL, NULL, NULL, NULL, NULL, NULL);
  dtrttp_(NULL, NULL, NULL, NULL, NULL, NULL);
  idamax_(NULL, NULL, NULL);
}