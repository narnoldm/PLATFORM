#include "metadata.hpp"
#include "param.hpp"

#include <set>

using namespace ::std;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    double t2, t1;

    if (rank != 0)
    {
        cout.rdbuf(sink.rdbuf());
    }

	// dimensions
	int M = 500;
	int N = 25;

	// set up matrices
	PGrid *evenG;
    evenG = new PGrid(rank, size, 0);
	pMat *A = new pMat(M, N, evenG, 0, 0, 0.0, false);
	pMat *b = new pMat(M, 1, evenG, 0, 0, 0.0, false);

	// fill matrices with random numbers
	srand(1);
	for (int i = 0; i < A->dataD.size(); ++i)
		A->dataD[i] = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);
	for (int i = 0; i < b->dataD.size(); ++i)
		b->dataD[i] = static_cast<double> (rand()) / static_cast<double> (RAND_MAX);

	// save to disk for checking
	A->write_bin("./A.bin");
	b->write_bin("./b.bin");

	// solve least squares
	b->leastSquares('N', M, N, 1, A, 0, 0, 0, 0);

	// write to disk
	b->write_bin("./x.bin");

	cout.rdbuf(strm_buffer);
    MPI_Finalize();
	

    return 0;
}