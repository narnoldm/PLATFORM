#ifndef SAMPLING_H
#define SAMPLING_H

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <mpi.h>
#include <set>

#include "pMat.hpp"
#include "processGrid.hpp"
#include "misc_func.hpp"

void random_oversampling(int nCells, int PointsNeeded, std::set<int>& samplingPoints);
void eigenvector_oversampling(pMat* URes, pMat* USol, int sampMethod, int nCells, int nVars, int nDOF, int numModesRHS,
							  int PointsNeeded, std::set<int>& samplingPoints, std::vector<int>& gP, std::string& timingOutput);
void gnat_oversampling_peherstorfer(pMat* URes, pMat* USol, int sampMethod, int nCells, int nVars, int nDOF, int numModesRHS,
							  int PointsNeeded, std::set<int>& samplingPoints, std::vector<int>& gP, std::string& timingOutput);
void gnat_oversampling_carlberg(pMat* URes, pMat* USol, int sampMethod, int nCells, int nVars, int nDOF, int numModesRHS,
							  int PointsNeeded, std::set<int>& samplingPoints, std::vector<int>& gP, std::string& timingOutput);

#endif