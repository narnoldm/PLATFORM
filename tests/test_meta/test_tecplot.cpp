

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

    cout << "testing tecplot metadata" << endl;

    

    string prefix = "testtec/gemsma_cmb_";
    string suffix = ".szplt";

    PGrid *loadG, *evenG;
    loadG = new PGrid(rank,size,1);
    evenG = new PGrid(rank,size,0);

    pMat *loadMat, *evenMatFromLoad, *evenMat;


    tecIO *dataset1;
    vector<string> token;
    token.push_back("input");
    token.push_back("tecplot");
    token.push_back(prefix);
    token.push_back(suffix);
    token.push_back("150000");
    token.push_back("151000");
    token.push_back("10");
    token.push_back("Static_Pressure");
    token.push_back("-1");
    token.push_back("Temperature");
    token.push_back("-1");
    dataset1 = new tecIO(token);
    string outdir="out";
    string outfile="stuff";


    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,loadG,0,0,0.0);
    dataset1->batchRead(loadMat);



    dataset1->batchWrite(loadMat,outdir,outfile);

    delete loadMat,evenMatFromLoad,evenMat;
    delete loadG,evenG,dataset1;


    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}