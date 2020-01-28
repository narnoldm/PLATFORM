
#ifndef EXTERN_FUNC_H
#define EXTERN_FUNC_H


extern "C" void Cblacs_pinfo(int *MYPNUM, int *NPROCS);
extern "C" void Cblacs_get(int ICONTXT, int WHAT, int *VAL);
extern "C" void Cblacs_gridinit(int *ICONTXT, char *ORDER, int prow, int pcol);
extern "C" void Cblacs_gridinfo(int ICONTXT, int *nprow, int *npcol, int *myprow, int *mypcol);
extern "C" void Cblacs_gridexit(int ICONTXT);
extern "C" double cblas_dnrm2(const MKL_INT n, const double *x, const MKL_INT incx);
extern "C" MKL_INT numroc_(const MKL_INT *n, const MKL_INT *nb, const MKL_INT *iproc, const MKL_INT *srcproc, const MKL_INT *nprocs);
extern "C" void descinit_(MKL_INT *desc, const MKL_INT *m, const MKL_INT *n, const MKL_INT *mb, const MKL_INT *nb, const MKL_INT *irsrc, const MKL_INT *icsrc, const MKL_INT *ictxt, const MKL_INT *lld, MKL_INT *info);
extern "C" void pdgeadd(const char *trans, const MKL_INT *m, const MKL_INT *n, const double *alpha, const double *a, const MKL_INT *ia, const MKL_INT *ja, const MKL_INT *desca, const double *beta, double *c, const MKL_INT *ic, const MKL_INT *jc, const MKL_INT *descc);
extern "C" void pdgemm(const char *transa, const char *transb, const MKL_INT *m, const MKL_INT *n, const MKL_INT *k, const double *alpha, const double *a, const MKL_INT *ia, const MKL_INT *ja, const MKL_INT *desca, const double *b, const MKL_INT *ib, const MKL_INT *jb, const MKL_INT *descb, const double *beta, double *c, const MKL_INT *ic, const MKL_INT *jc, const MKL_INT *descc);
extern "C" void pzgemm(const char *transa, const char *transb, const MKL_INT *m, const MKL_INT *n, const MKL_INT *k, const MKL_Complex16 *alpha, const MKL_Complex16 *a, const MKL_INT *ia, const MKL_INT *ja, const MKL_INT *desca, const MKL_Complex16 *b, const MKL_INT *ib, const MKL_INT *jb, const MKL_INT *descb, const MKL_Complex16 *beta, MKL_Complex16 *c, const MKL_INT *ic, const MKL_INT *jc, const MKL_INT *descc);
extern "C" void pdamax(const MKL_INT *n, double *amax, MKL_INT *indx, const double *x, const MKL_INT *ix, const MKL_INT *jx, const MKL_INT *descx, const MKL_INT *incx);
extern "C" void pdasum(const MKL_INT *n, double *asum, const double *x, const MKL_INT *ix, const MKL_INT *jx, const MKL_INT *descx, const MKL_INT *incx);
extern "C" void pdtran(const MKL_INT *m, const MKL_INT *n, const double *alpha, const double *a, const MKL_INT *ia, const MKL_INT *ja, const MKL_INT *desca, const double *beta, double *c, const MKL_INT *ic, const MKL_INT *jc, const MKL_INT *descc);


#endif