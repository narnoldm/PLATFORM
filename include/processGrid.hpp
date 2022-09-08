
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
	/// global number of process rows
	int prow;
	/// global number of process columns
	int pcol;
	/// dimensions of process grid
	int pdims[2];
	/// rank to be allowed to std out
	bool printRank;
	/// mpi rank 
	int rank;
	/// mpi size
	int size;

	PGrid(int, int, int); ///<Constructor of PBLACS Process Grid (int rank, int size, int type=0(rectangular),1(linear),2(1x))
	~PGrid(); ///<Destructor of PBLACS Process Grid
	int getDim(int); ///<Get the dimension of the process grid in the specified dimension
};

#endif