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
#include <assert.h>

//#include "extern_func.hpp"
#include "processGrid.hpp"



/*** 
*pMat: This file contains the headers for the PGrid and pMat classes.
*These define how the code distributes data into the ScaLAPACK format
* */


/***
 * pMat class is the underlying structure for all of PDP
 * contains the abstraction needed to call PBLACS and ScaLAPACK functions
 * */
class pMat
{
public:
	int myRC[2]={0,0}, desc[9]={0,0,0,0,0,0,0,0,0};
	int M=0, N=0, mb=0, nb=0;
	long MBs=0;
	std::vector<double> dataD;
#ifdef USE_MKL
	std::vector<MKL_Complex16> dataC;
#else
	std::vector<complex16> dataC;
#endif
	bool printRank=false;
	bool isComp=false;
	const int i_zero = 0;
	const int i_one = 1;

	int info=0, type=0, block=0, cycles=0;
	long long nelements=0;
	PGrid *pG;

	pMat();
	/// Will create copy pMat object of pointed one
	pMat(pMat *);
	pMat(int, int, PGrid *);
	pMat(int, int, PGrid *, bool);
	/// Creates pMat of dimension M,N on context pG with contant value
	pMat(int, int, PGrid *, int, int, double);
	pMat(int, int, PGrid *, int, int, double, bool);
	/// Creates pMat of dimension M,N on context pG with contant value with cycles
	pMat(int, int, PGrid *, int, int, int, double);
	pMat(int, int, PGrid *, int, int, int, double, bool);
	~pMat();

	/// core setup routine called by different contructors
	void setupMat(int, int, int, int, int, double, bool);
	/// Will swtich type
	void switchType(int);

	/// Will print out entire matrix (DO NOT USE unless debugging)
	void printMat();
	double getElement(int, int);
	double getLocalElement(int I, int J);

	//I/O
	int write_bin(std::string);
	int read_bin(std::string);
	bool check_bin_size(std::string, int &, int &);

	//PBLAS
	int matrix_Product(char tA, char tB, int m, int n, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc);
	int matrix_Sum(char, int, int, pMat *, int, int, int, int, double, double);
	int matrix_Product_sym(char uplo, char trans, int n, int k, double alpha, pMat *A, int ia, int ja, double beta, int ic, int jc);
	int matrix_vec_product(char trans, int m, int n, double alpha, pMat *A, int ia, int ja, pMat *B, int ib, int jb,
						   double beta, int ic, int jc);

	//Scalapack

	int svd_run(int, int, int, int, pMat *&, pMat *&, std::vector<double> &);
	int svd_run(int, int, int, int, pMat *&, pMat *&, std::vector<double> &, bool);

	//Utilities
	int transpose(pMat *);
	int transpose(pMat *, int, int, int, int);
	int changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb, bool stdout);
	int changeContext(pMat *, int, int, int, int, int, int);
	int changeContext(pMat *);
	int changeContext(pMat *, bool);
	int dMax(int, int, double &, int &);
	int dSum(int, int, double &);

	//Other
	int outerProductSum(pMat *U, char, pMat *VT, char, std::vector<double> &S, int inv);
	void pinv(pMat *A);
	int commCreate(MPI_Comm &col_com, int dim);
	
};

void destroyPMat(pMat *);
void destroyPMat(pMat *, bool);

std::ostream &operator<<(std::ostream &, const pMat &);
bool operator==(pMat const &, pMat const &);

#endif
