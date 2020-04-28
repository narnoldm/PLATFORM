
#ifndef EXTERN_FUNC_H
#define EXTERN_FUNC_H


/***
 * Ideally all external functions will be interfaced here
 * */

#include <mpi.h>


#ifdef USE_MKL
#include "mkl.h"
#include <mkl_cblas.h>
#include <mkl_pblas.h>
#include <mkl_scalapack.h>
#endif
//#include <pblas.h>







/// CBLACS functions
extern "C" void Cblacs_pinfo(int *MYPNUM, int *NPROCS);
extern "C" void Cblacs_get(int ICONTXT, int WHAT, int *VAL);
extern "C" void Cblacs_gridinit(int *ICONTXT, char *ORDER, int prow, int pcol);
extern "C" void Cblacs_gridinfo(int ICONTXT, int *nprow, int *npcol, int *myprow, int *mypcol);
extern "C" void Cblacs_gridexit(int ICONTXT);

#ifndef USE_MKL
/// BLAS routines


struct complex16
{
    double real;
    double imag;
};

extern "C" double cblas_dnrm2(const int n, const double *x, const int incx);

extern "C" int numroc_(const int *n, const int *nb, const int *iproc, const int *srcproc, const int *nprocs);
extern "C" void descinit_(int *desc, const int *m, const int *n, const int *mb, const int *nb, const int *irsrc, const int *icsrc, const int *ictxt, const int *lld, int *info);
extern "C" void pdgeadd(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pzgeadd(const char *trans, const int *m, const int *n, const complex16 *alpha, const complex16 *a, const int *ia, const int *ja, const int *desca, const complex16 *beta, complex16 *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *b, const int *ib, const int *jb, const int *descb, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdamax(const int *n, double *amax, int *indx, const double *x, const int *ix, const int *jx, const int *descx, const int *incx);
extern "C" void pdasum(const int *n, double *asum, const double *x, const int *ix, const int *jx, const int *descx, const int *incx);
extern "C" void pdtran(const int *m, const int *n, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pzgemm(const char *transa, const char *transb, const int *m, const int *n, const int *k, const complex16 *alpha, const complex16 *a, const int *ia, const int *ja, const int *desca, const complex16 *b, const int *ib, const int *jb, const int *descb, const complex16 *beta, complex16 *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdsyrk(const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdgemv(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *x, const int *ix, const int *jx, const int *descx, const int *incx, const double *beta, double *y, const int *iy, const int *jy, const int *descy, const int *incy);

extern "C" void pdgesvd(const char *,const char *, const int *,const int *,const double *,const int *,const int *,const int *, const double *, const double *,const int *,const int *,const int *,const double *,const int *,const int *,const int *, const double *,const int *,const int *);
extern "C" void pdgemr2d(const int *,const int *,const double *,const int *,const int *,const int *,const double *,const int *,const int *,const int *,const int *);
extern "C" void pzgemr2d(const int *,const int *,const complex16 *,const int *,const int *,const int *,const complex16 *,const int *,const int *,const int *,const int *);


#endif

#endif