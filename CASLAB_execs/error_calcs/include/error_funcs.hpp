#ifndef ERROR_FUNCS_H
#define ERROR_FUNCS_H

#include <string>
#include <vector>
#include <iostream>
#include <mpi.h>

#include "pMat.hpp"
#include "metadata.hpp"

void calc_integrated_error(pMat* dataMat, std::vector<std::string>& varNames, std::string outFile);
void calc_abs_and_l2_error(pMat* dataTruth, pMat* dataComp, tecIO* setData, std::string errDir, std::string errSuffix, bool outErrField);
void calc_abs_and_l2_error(pMat* dataTruth, pMat* dataComp, tecIO* setData, std::string errDir, std::string errSuffix, bool outErrField, bool mags, std::vector<std::string> magsToken);

#endif