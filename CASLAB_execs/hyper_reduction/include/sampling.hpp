#ifndef SAMPLING_H
#define SAMPLING_H

#include <string>
#include <vector>
#include <iostream>
#include <algorithm>
#include <mpi.h>
#include <unordered_set>

#include "pMat.hpp"
#include "param.hpp"
#include "processGrid.hpp"
#include "misc_func.hpp"

void qr_sampling(paramMap inputFile, const std::string& qrSampFileStr, const std::string& outFileStr, int nCells,
                 pMat* U_T, std::vector<int>& gP, std::unordered_set<int>& samplingPoints);
void random_oversampling(int nCells, int PointsNeeded, std::unordered_set<int>& samplingPoints, std::vector<int>& gP);

void eigenvector_oversampling(std::vector<pMat*> U_vec, int sampMethod, int nCells, int nVars, int PointsNeeded,
                              std::unordered_set<int>& samplingPoints, std::vector<int>& gP, std::string& timingOutput);
void eigenvector_oversampling_metric(pMat* U, pMat* U_samp, pMat* U_samp_copy, pMat* rVec, pMat* rVecSum, pMat* nonUniqueVec, int numCurrentDOFs);

void gnat_oversampling_peherstorfer(std::vector<pMat*> U_vec, int sampMethod, int nCells, int nVars, int PointsNeeded,
									std::unordered_set<int>& samplingPoints, std::vector<int>& gP, std::string& timingOutput);
void gnat_oversampling_peherstorfer_metric(pMat* U, pMat* U_samp, pMat* lsSol, int iterNum, pMat* rVec, pMat* rVecSum, pMat* nonUniqueVec, int numCurrentDOFs);

void gnat_oversampling_carlberg(std::vector<pMat*> U_vec, int sampMethod, int nCells, int nVars, int PointsNeeded,
								std::unordered_set<int>& samplingPoints, std::vector<int>& gP, std::string& timingOutput);
void gnat_oversampling_carlberg_metric(pMat* U, pMat* U_samp, pMat* lsSol, int iterNum, pMat* rVec, pMat* rVecSum, pMat* nonUniqueVec, int numCurrentDOFs);

void calc_regressor(pMat* U, pMat* regressor, std::vector<int>& gP, int nCells, int nVars);
void emplace_zeros(const char transA, pMat* AIn, pMat* AOut, std::vector<int>& gP, int nCells, int nVars);

#endif
