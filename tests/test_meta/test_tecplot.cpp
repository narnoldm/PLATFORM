

#include "metadata.hpp"

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

    cout << "testing metadata" << endl;

    meta *dataset1;

    string prefix = "testplit/test";
    string suffix = ".bin";

    dataset1 = new meta(1,100,1,prefix,suffix);

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}