#include "metadata.hpp"
#include "param.hpp"

using namespace ::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();


    int debug_proc = 0;

    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }


    string input = "";
    input = argv[1];

    cout << "FileName is: " << input << endl;

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);

    pMat *evenMat;

    int nPoints=0, nSets=0;
    evenMat->check_bin_size(input, nPoints, nSets);
    evenMat = new pMat(nPoints, nSets, evenG, 0, 0, 0.0);

    evenMat->read_bin(input);

    pMat *V, *YT;
    vector<double> S;

    int minDim = min(nPoints, nSets);

    V = new pMat(nPoints, minDim, evenG, 0, 0, 0.0);
    YT = new pMat(minDim, nSets, evenG, 0, 0, 0.0);

    S.resize(minDim);

    evenMat->svd_run(nPoints, nSets, 0, 0, V, YT, S);

    printASCIIVecP0("S.txt", S, S.size());

    delete evenMat;

    YT->write_bin("YT.bin");
    V->write_bin("V.bin");

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}