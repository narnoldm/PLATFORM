#ifndef MISC_FUNC_MPI_H
#define MISC_FUNC_MPI_H

#include <string>
#include <fstream>
#include <cmath>
#include <mpi.h>

void aggregateTiming(double t, std::string filename, std::string label);

#endif