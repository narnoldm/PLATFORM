#include "misc_func_mpi.hpp"

using namespace :: std;

/**
 * Computes aggregate timing statistics across all MPI processes. 
 * 
 * Calculates average, minimum, maximum, and standard deviation.
 * Appends recorded results to end of log file given by filename, with label given by label
*/
void aggregateTiming(double t, string filename, string label) {

	int rank, num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// average
	double average;
	MPI_Allreduce(&t, &average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	average /= num_procs;
     
    // min and max
	double maximum, minimum;
    MPI_Reduce(&t, &maximum, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t, &minimum, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD);

    // standard deviation
	double global_squareDiff;
    double local_squareDiff = pow(t - average, 2.0);
    MPI_Reduce(&local_squareDiff, &global_squareDiff, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double standard_dev = sqrt(global_squareDiff / double(num_procs - 1));
		ofstream out;
		out.open(filename, ios::app);
		out << "=================================" << endl;
		out << label << " timings" << endl;
		out << "=================================" << endl;
		out << "Average: " << to_string(average) << endl;
		out << "Maximum: " << to_string(maximum) << endl;
		out << "Minimum: " << to_string(minimum) << endl;
		out << "STD: " << to_string(standard_dev) << endl;
		out.close();
	}

}
