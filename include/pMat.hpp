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

#include "extern_func.hpp"

/*** 
*pMat: This file contains the headers for the PGrid and pMat classes.
*These define how the code distributes data into the ScaLAPACK format
* */

/***
 * PGrid classes are the containers that have the PBLACS contexts
 * and communicators. All PBLAS distributed matrixes need this information
 * All pMats are built on an underlying Process Grid (PGrid). Multiple 
 *	PMats can be assigned to identical PGrids. 
 * */

class PGrid
{
public:
	/// Context Identifier
	int icntxt;
	/// local process row index
	int myrow;
	/// local process col identifier
	int mycol;
	/// global number of rows
	int prow;
	/// global number of columns
	int pcol;
	/// dimesnions of process grid
	int pdims[2];
	/// rank to be allowed to std out
	bool printRank;
	/// mpi rank and size
	int rank, size;

	PGrid(int r, int s, int type); //Constructor
	~PGrid();
	int getDim(int dim);
};
/***
 * pMat class is the underlying structure for all of PDP
 * contains the abstraction needed to call PBLACS and ScaLAPACK functions
 * Refactor is currently needed (note N is number of rows M is number of columns)
 * */
class pMat
{
public:
	int myRC[2], desc[9];
	int M, N, mb, nb;
	long MBs;
	std::vector<double> dataD;
	#ifdef USE_MKL
	std::vector<MKL_Complex16> dataC;
	#else
	std::vector<complex16> dataC;
	#endif
	bool printRank;
	bool isComp;
	const int i_zero = 0;
	const int i_one = 1;

	int info, type, block, cycles;
	long long nelements;
	PGrid *pG;

	pMat();
	/// Will create copy pMat object of pointed one
	pMat(pMat *point);
	/// Creates pMat of dimension M,N on context pG with contant value
	pMat(int m, int n, PGrid *pG, int t, int b, double init);
	/// Creates pMat of dimension M,N on context pG with contant value with cycles
	pMat(int m, int n, PGrid *pG, int t, int b, int c, double init);
	~pMat();

	/// core setup routine called by different contructors
	void setupMat(int m, int n, int t, int b, int c, double init);
	/// Will swtich type
	void switchType(int t);

	/// Will print out entire matrix (DO NOT USE unless debugging)
	void printMat();

	//I/O 
	int write_bin(std::string filename);
	int read_bin(std::string &filename);
	bool check_bin_size(std::string filename, int &mN, int &mM);

	//PBLAS
	int matrix_Product(char tA, char tB, int m, int n, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc);
	int matrix_Sum(char tA, int m, int n, pMat *A, int ia, int ja, int ib, int jb, double alpha, double beta);

	//Scalapack

	int svd_run(int m, int n, int ia, int ja, pMat *&U, pMat *&VT, std::vector<double> &S);

	//Utilities
	int transpose(pMat *A, int m, int n, int ia, int ja);
	int changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb);
	int changeContext(pMat *A);
	int dMax(int dim, int rc, double &val);
	int dAve(int dim, int rc, double &val);
};

std::ostream &operator<<(std::ostream &os, const pMat &p);
bool operator==(pMat const &p1, pMat const &p2);

#endif
