

#include "metadata.hpp"



using namespace::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();
    int out=0;
    if(argc==2)
        out=atoi(argv[1]);
    
    if (rank != out)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    PGrid *evenG;
    evenG = new PGrid(rank,size,0);

    

    string prefix = "testsplit/test";
    string suffix = ".bin";


    pMat *loadMat, *evenMatFromLoad, *evenMat;


    meta *dataset1;
    dataset1 = new meta(1,11,1,prefix,suffix);
    string outdir="out";
    string outfile="stuff";


    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,1,0.0);
    dataset1->batchRead(loadMat);

    evenMat=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);

    evenMat->changeContext(loadMat);
    
    evenMatFromLoad = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);

    dataset1->batchRead(evenMatFromLoad);

    evenMatFromLoad->printMat();

    assert((*evenMat)==(*evenMatFromLoad));   



    dataset1->batchWrite(evenMatFromLoad,"out1","outtest");

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}