#ifndef PMAT_H
#define PMAT_H

#include <stdlib.h>

#include <vector>
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <cmath>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <limits>

#include "extern_func.hpp" // TODO: this needs a compiler flag to avoid adding when not using MKL
#include "processGrid.hpp"

/*#ifdef USE_MKL
#ifndef EIGEN_USE_MKL_ALL
#define EIGEN_USE_MKL_ALL
#endif
#endif*/
#include <Eigen/Dense>

#ifndef MAX_OPERATION_SIZE
#define MAX_OPERATION_SIZE 10000
#endif

/***
*pMat: This file contains the headers for the PGrid and pMat classes.
*These define how the code distributes data into the ScaLAPACK format
* */


/***
 * pMat class is the underlying structure for all of PDP
 * contains the abstraction needed to call PBLACS and ScaLAPACK functions
 * Refactor is currently needed (note N is number of rows M is number of columns)
 * */
class pMat
{
public:
	int myRC[2] = {0,0};
    int* desc;
	int M = 0, N = 0, mb = 0, nb = 0;
	long MBs = 0;

	std::vector<double> dataD;
#ifdef USE_MKL
	std::vector<MKL_Complex16> dataC;
#else
	std::vector<complex16> dataC;
#endif

	bool printRank = false;
	bool isComp = false;
	const int i_zero = 0;
	const int i_one = 1;

	int info = 0, type = 0, block = 0, cycles = 0;
	long long nelements = 0;
	PGrid *pG;

	pMat();
	// Will create copy pMat object of pointed one
	pMat(pMat *point);
	pMat(int m, int n, PGrid *pGp);
	pMat(int m, int n, PGrid *pGp, bool);
	// Creates pMat of dimension M,N on context pG with contant value
	pMat(int m, int n, PGrid *pGp, int t, int b, double);
	pMat(int m, int n, PGrid *pGp, int t, int b, double, bool stdout_bool);
	// Creates pMat of dimension M,N on context pG with contant value with cycles
	pMat(int m, int n, PGrid *pGp, int t, int b, int c, double);
	pMat(int m, int n, PGrid *pGp, int t, int b, int c, double, bool stdout_bool);
	~pMat();

	// core setup routine called by different contructors
	void setupMat(int m, int n, int t, int b, int c, double init, bool stdout);

    // Will switch type
	void switchType(int t);

	// Will print out entire matrix (DO NOT USE unless debugging)
	void printMat();

	// I/O
	void write_bin(std::string filename);
	void read_bin(std::string filename);
	bool check_bin_size(std::string filename, int &mM, int &mN);

	// PBLAS
	void matrix_Product(char tA, char tB, int m, int n, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc);
	void matrix_Sum(char tA, int m, int n, pMat *A, int ia, int ja, int ib, int jb, double alpha, double beta);
    void ATA_lowMem(int m, int n, double alpha, pMat* A, int ia, int ja, int ic, int jc);
	void matrix_Product_sym(char uplo, char trans, int n, int k, double alpha, pMat *A, int ia, int ja, double beta, int ic, int jc);
	void matrix_vec_product(char trans, int m, int n, double alpha, pMat *A, int ia, int ja, pMat *B, int ib, int jb,
						   double beta, int ic, int jc);
    void scale_col_row(double alpha, int idx, bool scaleRows);
    void matrix_elem_mult(char tA, int m, int n, double alpha, pMat *A, int ia, int ja, int ic, int jc);

	// SVD
	void svd_run(int, int, int, int, pMat *&, pMat *&, std::vector<double> &);
	void svd_run(int, int, int, int, pMat *&, pMat *&, std::vector<double> &, bool);

	// MOS
	void mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, std::vector<double> &S);
	void mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, std::vector<double> &S, int modeStart, int modeEnd);
	void mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, std::vector<double> &S, int modeStart, int modeEnd,
				int mosStep, PGrid *procGrid);

	// QR decomposition
	void qr_run(int n, int m, int ia, int ja, std::vector<int> &ipiv);
	void qr_run(int n, int m, int ia, int ja, std::vector<int> &ipiv, std::string outdir, bool stdout);
	void qr_run(int n, int m, int ia, int ja, std::vector<int> &ipiv, std::string outdir, std::string outfile, bool stdout);

	// Least-squares
	void leastSquares(char trans, int m, int n, int nrhs, pMat *&A, int ia, int ja, int ib, int jb);

	// Utilities
	void transpose(pMat *);
	void transpose(pMat *, int, int, int, int);
	void changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb, bool stdout);
	void changeContext(pMat *, int, int, int, int, int, int);
	void changeContext(pMat *);
	void changeContext(pMat *, bool);
	void dMax(int, int, double &, int &);
	int argmax_vec();
	void dSum(int, int, double &);
    double getElement(int I, int J);
	double getLocalElement(int I, int J);
    double getLocalElement(int I, int J, double temp);
	void setElement(int I, int J, double val);
    void setElement(int I, int J, double val, bool barrier);
    int getDataIndex(int I, int J);

	// Other
	void outerProductSum(pMat *U, char, pMat *VT, char, std::vector<double> &S, int inv);
	void pinv(pMat *A);
	void commCreate(MPI_Comm &col_com, int dim);

};

void destroyPMat(pMat *);
void destroyPMat(pMat *, bool);

std::ostream &operator<<(std::ostream &, const pMat &);
bool operator==(pMat const &, pMat const &);

#endif
