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


#include "extern_func.hpp"


using namespace ::std;





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
	
	//I/O
	int write_bin(std::string filename);
	int read_bin(string &filename);
	bool check_bin_size(string filename, int &mN, int &mM);
	
	//PBLAS
	int matrix_Product(char tA, char tB, int m, int n, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc);
	int matrix_Sum(char tA, int m, int n, pMat *A, int ia, int ja, int ib, int jb, double alpha, double beta);
	
	//Scalapack
	
	int svd_run(int m, int n, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S);
	
	//Utilities
	int transpose(pMat *A, int m, int n, int ia, int ja);
	int changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb);
	int changeContext(pMat *A);
	int dMax(int dim, int rc, double &val);
	int dAve(int dim, int rc, double &val);
};

ostream &operator<<(std::ostream &os, const pMat &p);
bool operator==(pMat const &p1, pMat const &p2);
