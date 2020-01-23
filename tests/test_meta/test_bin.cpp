

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

    

    string prefix = "testsplit/test";
    string suffix = ".bin";

    PGrid *loadG, *evenG;
    loadG = new PGrid(rank,size,1);
    evenG = new PGrid(rank,size,0);

    pMat *loadMat, *evenMatFromLoad, *evenMat;


    meta *dataset1;
    dataset1 = new meta(1,100,1,prefix,suffix);

    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,loadG,0,0,0.0);
    dataset1->batchRead(loadMat);

    evenMatFromLoad = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    evenMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);

    string evenMatFN= "test.bin";
    evenMat->read_bin(evenMatFN);


    evenMatFromLoad->changeContext(loadMat);


    assert((*evenMatFromLoad)==(*evenMat));   




    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}