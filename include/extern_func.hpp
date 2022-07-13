
#ifndef EXTERN_FUNC_H
#define EXTERN_FUNC_H

/***
 * Ideally all external functions will be interfaced here
 * */

#include <mpi.h>
#include <string.h>
#include <algorithm>
#include <iostream>

#ifdef USE_MKL
#include "mkl.h"
#include <mkl_cblas.h>
#include <mkl_pblas.h>
#include <mkl_scalapack.h>

#endif
//#include <pblas.h>

// CBLACS functions
extern "C" void Cblacs_pinfo(int *MYPNUM, int *NPROCS);
extern "C" void Cblacs_get(int ICONTXT, int WHAT, int *VAL);
extern "C" void Cblacs_gridinit(int *ICONTXT, char *ORDER, int prow, int pcol);
extern "C" void Cblacs_gridinfo(int ICONTXT, int *nprow, int *npcol, int *myprow, int *mypcol);
extern "C" void Cblacs_gridexit(int ICONTXT);

#ifndef USE_MKL

// PTOOLS routines
extern "C"
{
    struct PB_VM_T
    {
        int offd;
        int lcmt00;
        int mp;
        int imb1;
        int imbloc;
        int mb;
        int lmbloc;
        int mblks;
        int iupp;
        int upp;
        int prow;
        int nprow;
        int nq;
        int inb1;
        int inbloc;
        int nb;
        int lnbloc;
        int nblks;
        int ilow;
        int low;
        int pcol;
        int npcol;
        int lcmb;
    };

    typedef char * F_CHAR_T;

    typedef void (*GESD2D_T)   (int, int, int, char *, int, int, int);
    typedef void (*GERV2D_T)   (int,      int,       int, char *,    int,       int, int );
    typedef void (*GEBS2D_T)   (int,      char *,    char *, int,       int,       char *,int );
    typedef void (*GEBR2D_T)   (int,      char *,    char *,int,       int,       char *,int,       int,       int );
    typedef void (*GSUM2D_T)   (int,      char *,    char *,int,       int,       char *,int,       int,       int );
    typedef void (*MMADD_T)    (int *, int *, char *, char *, int*, char*, char*, int*);
    typedef void (*MMSHFT_T)   (int  *,   int  *,    int *, char *,    int  * );
    typedef void (*VVDOT_T)    (int  *,   char *,    char *, int  *,    char *,    int  * );
    typedef void (*VVSET_T)    (int  *,   char *,    char *, int  * );
    typedef void (*TZPAD_T)    (F_CHAR_T, F_CHAR_T,  int  *, int  *,    int  *,    char *, char *,    char *,    int  * );
    typedef void (*TZPADCPY_T) (F_CHAR_T, F_CHAR_T,  int  *,int  *,    int  *,    char *,int *,     char *,    int  * );
    typedef void (*TZSET_T)    (F_CHAR_T, int  *,    int  *,int  *,    char *,    char *,char *,    int  * );
    typedef void (*TZSCAL_T)   (F_CHAR_T, int *,     int  *,int  *,    char *,    char *,int  * );
    typedef void (*AXPY_T)     (int *,    char *,    char *,int *,     char *,    int * );
    typedef void (*COPY_T)     (int *,    char *,    int *,char *,    int * );
    typedef void (*SWAP_T)     (int *,    char *,    int *, char *,    int * );
    typedef void (*GEMV_T)     (F_CHAR_T, int*, int*, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*AGEMV_T)    (F_CHAR_T, int*, int*, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*SYMV_T)     (F_CHAR_T, int*, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*ASYMV_T)    (F_CHAR_T, int*, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*HEMV_T)     (F_CHAR_T, int*, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*AHEMV_T)    (F_CHAR_T, int*, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*TRMV_T)     (F_CHAR_T, F_CHAR_T, F_CHAR_T, int*, char*, int*, char*, int*);
    typedef void (*ATRMV_T)    (F_CHAR_T, F_CHAR_T, F_CHAR_T, int*, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*TRSV_T)     (F_CHAR_T, F_CHAR_T, F_CHAR_T, int*, char*, int*, char*, int*);
    typedef void (*GERC_T)     (int*, int*, char*, char*, int*, char*, int*, char*, int*);
    typedef void (*GERU_T)     (int*, int*, char*, char*, int*, char*, int*, char*, int*);
    typedef void (*SYR_T)      (F_CHAR_T, int*, char*, char*, int*, char*, int*);
    typedef void (*HER_T)      (F_CHAR_T, int*, char*, char*, int*, char*, int*);
    typedef void (*SYR2_T)     (F_CHAR_T, int*, char*, char*, int*, char*, int*, char*, int*);
    typedef void (*HER2_T)     (F_CHAR_T, int*, char*, char*, int*, char*, int*, char*, int*);
    typedef void (*GEMM_T)     (F_CHAR_T, F_CHAR_T, int*, int*,int *, char*, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*SYMM_T)     (F_CHAR_T, F_CHAR_T, int*, int*,char *, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*HEMM_T)     (F_CHAR_T, F_CHAR_T, int*, int*,char *, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*SYRK_T)     (F_CHAR_T, F_CHAR_T, int*, int*,char *, char*, int*, char*, char*, int*);
    typedef void (*HERK_T)     (F_CHAR_T, F_CHAR_T, int*, int*,char *, char*, int*, char*, char*, int*);
    typedef void (*SYR2K_T)    (F_CHAR_T, F_CHAR_T, int*, int*,char *, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*HER2K_T)    (F_CHAR_T, F_CHAR_T, int*, int*,char *, char*, int*, char*, int*, char*, char*, int*);
    typedef void (*TRMM_T)     (F_CHAR_T, F_CHAR_T, F_CHAR_T, F_CHAR_T, int*, int*, char*, char*, int*, char*, int*);
    typedef void (*TRSM_T)     (F_CHAR_T, F_CHAR_T, F_CHAR_T, F_CHAR_T, int*, int*, char*, char*, int*, char*, int*);

    struct PBTYP_T
    {
        char type;
        int  usiz;
        int  size;

        char * zero, * one, * negone;

        GESD2D_T  Cgesd2d;
        GERV2D_T  Cgerv2d;
        GEBS2D_T  Cgebs2d;
        GEBR2D_T  Cgebr2d;
        GSUM2D_T  Cgsum2d;
        MMADD_T  Fmmadd;
        MMADD_T  Fmmcadd;
        MMADD_T  Fmmtadd;
        MMADD_T  Fmmtcadd;
        MMADD_T  Fmmdda;
        MMADD_T  Fmmddac;
        MMADD_T  Fmmddat;
        MMADD_T  Fmmddact;
        MMSHFT_T Fcshft;
        MMSHFT_T Frshft;
        VVDOT_T  Fvvdotu;
        VVDOT_T  Fvvdotc;
        TZPAD_T  Ftzpad;
        TZPADCPY_T Ftzpadcpy;
        VVSET_T Fset;
        TZSCAL_T Ftzscal;
        TZSCAL_T Fhescal;
        TZSCAL_T Ftzcnjg;
        AXPY_T Faxpy;
        COPY_T Fcopy;
        SWAP_T Fswap;
        GEMV_T Fgemv;
        SYMV_T Fsymv;
        HEMV_T Fhemv;
        TRMV_T Ftrmv;
        TRSV_T Ftrsv;
        AGEMV_T Fagemv;
        ASYMV_T Fasymv;
        AHEMV_T Fahemv;
        ATRMV_T Fatrmv;
        GERC_T Fgerc;
        GERU_T Fgeru;
        SYR_T Fsyr;
        HER_T Fher;
        SYR2_T Fsyr2;
        HER2_T Fher2;
        GEMM_T Fgemm;
        SYMM_T Fsymm;
        HEMM_T Fhemm;
        SYRK_T Fsyrk;
        HERK_T Fherk;
        SYR2K_T Fsyr2k;
        HER2K_T Fher2k;
        TRMM_T Ftrmm;
        TRSM_T Ftrsm;

    };

};

#define F2C_CHAR(a) (a)
#define C2F_CHAR(a) (a)
#define Mptr( a_, i_, j_, lda_, siz_ ) ( (a_) + ( (long long) ( (long long)(i_)+ (long long)(j_)*(long long)(lda_))*(long long)(siz_) ) )
#define MModAdd(I1, I2, d) ( ( (I1) + (I2) < (d) ) ? (I1) + (I2) : (I1) + (I2) - (d) )
#define MModAdd1(I, d) ( ((I) != (d)-1) ? (I) + 1 : 0 )
#define MModSub(I1, I2, d) ( ( (I1) < (I2) ) ? (d) + (I1) - (I2) : (I1) - (I2) )
#define MModSub1(I, d) ( ((I)!=0) ? (I)-1 : (d)-1 )
#define MPosMod(I, d) ( (I) - ((I)/(d))*(d) )
#define Mupcase(C) (((C)>96 && (C)<123) ? (C) & 0xDF : (C))
#define DTYPE_ 0
#define CTXT_  1
#define M_     2
#define N_     3
#define IMB_   4
#define INB_   5
#define MB_    6
#define NB_    7
#define RSRC_  8
#define CSRC_  9
#define LLD_   10
#define DLEN_  11
#define CCOLUMN    'C'
#define CROW       'R'
#define CNOTRAN    'N'
#define CNOCONJG   'N'
#define CTRAN      'T'
#define CCONJG     'Z'
#define CCOTRAN    'C'
#define CPACKING   'P'
#define CUNPACKING 'U'


extern "C" PBTYP_T* PB_Cdtypeset();
extern "C" int PB_Cspan(int N, int I, int INB, int NB, int SRCPROC, int NPROCS);
extern "C" void PB_CargFtoC(int IF, int JF, int *DESCIN, int *IC, int *JC, int *DESCOUT);
extern "C" int PB_Cnumroc(int N, int I, int INB, int NB, int PROC, int SRCPROC, int NPROCS);
extern "C" int PB_Cfirstnb(int N, int I, int INB, int NB);
extern "C" int pilaenv_(const int *ICTXT, const char *PREC);
extern "C" int PB_Clcm(int M, int N);
extern "C" int PB_Cgcd(int M, int N);
extern "C" int PB_Cindxg2p(int IG, int INB, int NB, int PROC, int SRCPROC, int NPROCS);
extern "C" int PB_CVMnpq(PB_VM_T *VM);
extern "C" char* PB_Cmalloc(int LENGTH);
extern "C" void PB_CVMcontig(PB_VM_T *VM, int *NRPQ, int *NCPQ, int *IOFF, int *JOFF);
extern "C" void PB_Cdescset(int *DESC, int M, int N, int IMB, int INB, int MB, int NB,
    int RSRC, int CSRC, int CTXT, int LLD);
extern "C" void PB_CVMupdate(PB_VM_T *VM, int K, int *II, int *JJ);
extern "C" void PB_CVMinit(PB_VM_T *VM, int OFFD, int M, int N, int IMB1, int INB1,
    int MB, int NB, int MRROW, int MRCOL, int NPROW, int NPCOL, int LCMB);
extern "C" void PB_Cinfog2l(int I, int J, int* DESC, int NPROW, int NPCOL,
    int MYROW, int MYCOL, int* II, int* JJ, int* PROW, int* PCOL);
extern "C" void PB_Cpaxpby(PBTYP_T * TYPE, char * CONJUG, int M, int N,
    char * ALPHA, char * A, int IA, int JA, int * DESCA, char * AROC,
    char * BETA, char * B, int IB, int JB, int * DESCB, char * BROC);
extern "C" int PB_CVMpack(PBTYP_T *TYPE, PB_VM_T *VM, char *VROCS, char *ROCS,
    char *UNPA, char *TRANS, int MN, int K, char *ALPHA, char *A, int LDA,
    char *BETA, char *B, int LDB);
extern "C" int PB_CVMloc(PBTYP_T* TYPE, PB_VM_T* VM, char* VROCS,
    char* ROCS, char* UNPA, char* TRANS, int MN, int K, char* ALPHA,
    char* A, int LDA, char* BETA,  char* B, int LDB);
extern "C" char * PB_Ctop(int* ICTXT, char* OP, char* SCOPE, char* TOP);

// custom functions
void pdgemul(char trans, int m, int n, double alpha,
    double* AMat, int ai, int aj, int* adesc,
    double* CMat, int ci, int cj, int* cdesc);
void daxy(PBTYP_T* TYPE, char* CONJUG, int M, int N, char* ALPHA,
    char* A, int IA, int JA, int* DESCA, char* AROC,
    char* B, int IB, int JB, int* DESCB, char* BROC);
int VMpackMult(PBTYP_T* TYPE, PB_VM_T* VM, char* VROCS, char* ROCS,
    char* UNPA, char* TRANS, int MN, int K, char* ALPHA,
    char* A, int LDA,
    char* B, int LDB);
int VMlocMult( PBTYP_T * TYPE, PB_VM_T * VM, char * VROCS, char * ROCS,
    char * UNPA, char * TRANS, int MN, int K, char * ALPHA,
    char * A, int LDA,
    char * B, int LDB );

typedef void (*MMMULT_T) (int, int, char*, char*, int, char*, int);

// BLAS routines
struct complex16
{
    double real;
    double imag;
};

extern "C" double cblas_dnrm2(const int n, const double *x, const int incx);

extern "C" int numroc_(const int *n, const int *nb, const int *iproc, const int *srcproc, const int *nprocs);
extern "C" void descinit_(int *desc, const int *m, const int *n, const int *mb, const int *nb, const int *irsrc, const int *icsrc, const int *ictxt, const int *lld, int *info);
extern "C" void pdgeadd_(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pzgeadd_(const char *trans, const int *m, const int *n, const complex16 *alpha, const complex16 *a, const int *ia, const int *ja, const int *desca, const complex16 *beta, complex16 *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *b, const int *ib, const int *jb, const int *descb, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdamax_(const int *n, double *amax, int *indx, const double *x, const int *ix, const int *jx, const int *descx, const int *incx);
extern "C" void pdasum_(const int *n, double *asum, const double *x, const int *ix, const int *jx, const int *descx, const int *incx);
extern "C" void pdtran_(const int *m, const int *n, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pzgemm_(const char *transa, const char *transb, const int *m, const int *n, const int *k, const complex16 *alpha, const complex16 *a, const int *ia, const int *ja, const int *desca, const complex16 *b, const int *ib, const int *jb, const int *descb, const complex16 *beta, complex16 *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdsyrk_(const char *uplo, const char *trans, const int *n, const int *k, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *beta, double *c, const int *ic, const int *jc, const int *descc);
extern "C" void pdgemv_(const char *trans, const int *m, const int *n, const double *alpha, const double *a, const int *ia, const int *ja, const int *desca, const double *x, const int *ix, const int *jx, const int *descx, const int *incx, const double *beta, double *y, const int *iy, const int *jy, const int *descy, const int *incy);
extern "C" void pdscal_(const int *n, const double *alpha, const double *x, const int *ix, const int *jx, const int *desc, const int *incx);

extern "C" void pdgesvd_(const char *, const char *, const int *, const int *, const double *, const int *, const int *, const int *, const double *, const double *, const int *, const int *, const int *, const double *, const int *, const int *, const int *, const double *, const int *, const int *);
extern "C" void pdgesvd_(const char *, const char *, const int *, const int *, const double *, const int *, const int *, const int *, const double *, const double *, const int *, const int *, const int *, const double *, const int *, const int *, const int *, const double *, const int *, const int *);
extern "C" void pdgeqpf_(const int *, const int *, double *, const int *, const int *, const int *, const int *, double *, double *work, const int *, const int *);
extern "C" void pdgels_(const char *trans, const int *m, const int *n, const int *nrhs, const double *a, const int *ia, const int *ja, const int *desca, const double *b, const int *ib, const int *jb, const int *descb, const double *work, const int *lwork, const int *info);

extern "C" void pdgemr2d_(const int *, const int *, const double *, const int *, const int *, const int *, const double *, const int *, const int *, const int *, const int *);
extern "C" void pzgemr2d_(const int *, const int *, const complex16 *, const int *, const int *, const int *, const complex16 *, const int *, const int *, const int *, const int *);

#endif

#endif
