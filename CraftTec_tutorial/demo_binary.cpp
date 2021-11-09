

#include "metadata.hpp"
#include "param.hpp"

using namespace::std;
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

    paramMap p1("demoInput.txt",rank);

    string prefix = "testsplit/test";
    string suffix = ".bin";

    PGrid *evenG;
    evenG = new PGrid(rank,size,0);

    pMat *evenMat;

    string outdir="out";
    string outfile="U";
    p1.getParamString("outdir",outdir);
    p1.getParamString("outfile",outfile);

    meta *dataset1;
    vector<string> token;

    string lToken;
    p1.getParamString("token",lToken);

 
    tokenparse(lToken,"|",token);
    dataset1 = new meta(token);


    evenMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    dataset1->batchRead(evenMat);


    pMat *U,*VT;
    vector<double> S;

    U=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    VT=new pMat(dataset1->nSets,dataset1->nSets,evenG,0,0,0.0);
    S.resize(dataset1->nSets);
    evenMat->svd_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S);
    evenMat->outerProductSum(U,'N',VT,'N',S,0);


    dataset1->batchWrite(U,outdir,outfile);

   


    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}