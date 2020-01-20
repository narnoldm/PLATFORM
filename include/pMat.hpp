#include <stdlib.h>
#include <mpi.h>
#include "mkl.h"
#include <mkl_cblas.h>
#include <mkl_pblas.h>
#include <mkl_scalapack.h>
#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <sstream>

using namespace ::std;

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

/* 
pMat: This file contains the headers for the PGrid and pMat classes.
These define how the code distributes data into the ScaLAPACK format
 */

class PGrid
{
/* 
	PGrid sets up the communication context needed by all PBLACS routines.
	All pMats are built on an underlying Process Grid (PGrid). Multiple 
	PMats can be assigned to identical PGrids. 
 */
public:
	int icntxt; //context identifier
	int myrow; //(Local)
	int mycol; 
	int prow;
	int pcol;
	int pdims[2];
	bool printRank; 
	int rank, size;

	PGrid(int r, int s, int type); //Constructor
	~PGrid();
	int getDim(int dim);
};

class pMat
{
public:
	int myRC[2], desc[9];
	int M, N, mb, nb;
	long MBs;
	std::vector<double> dataD;
	std::vector<MKL_Complex16> dataC;
	bool printRank;
	bool isComp;
	const int i_zero = 0;
	const int i_one = 1;

	int info, type, block, cycles;
	long long nelements;
	PGrid *pG;

	pMat();
	pMat(pMat * point);
	pMat(int m, int n, PGrid *pG, int t, int b, double init);
	pMat(int m, int n, PGrid *pG, int t, int b, int c, double init);
	~pMat();
	void setupMat(int m, int n, int t, int b, int c, double init);
	void switchType(int t);
	void printMat();
	int write_bin(std::string filename);
	int read_bins(std::string prefix, int start, int skip, int end);
	int read_single_bin(const char *name, int col);
	int read_bin(string &filename);
	bool check_bin_size(string filename, int &mN, int &mM);
	int matrix_Product(char tA, char tB, int m, int n, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc);
	int matrix_Sum(char tA, int m, int n, pMat *A, int ia, int ja, int ib, int jb, double alpha, double beta);
	int svd_run(int m, int n, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S);
	int transpose(pMat *A, int m, int n, int ia, int ja);
	int changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb);
	int changeContext(pMat *A);
	int dMax(int dim, int rc, double val);
	int dAve(int dim, int rc, double val);
	void printDesc();
};

ostream &operator<<(std::ostream &os, const pMat &p);
bool operator==(pMat const &p1, pMat const &p2);
