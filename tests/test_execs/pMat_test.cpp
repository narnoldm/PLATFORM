

#include "pMat.hpp"
#include <assert.h>


int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();
    if (rank != 0)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    PGrid *p1 = new PGrid(rank,size,0);
    pMat *m1;
    int M,N;
    string filename="test.bin";
    assert(m1->check_bin_size(filename,M,N));
    m1=new pMat(M,N,p1,1,1,0.0);
    m1->read_bin(filename);




    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}