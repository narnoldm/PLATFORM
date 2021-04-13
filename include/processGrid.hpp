
#ifndef PROCESSGRID_H
#define PROCESSGRID_H

#include "extern_func.hpp"
#include <iostream>
#include <string>



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

	PGrid(int, int, int); //Constructor
    PGrid(int, int);
	~PGrid();
	int getDim(int);
};

#endif