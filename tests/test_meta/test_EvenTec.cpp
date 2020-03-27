

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
   

    string prefix = "testtec/gemsma_cmb_";
    string suffix = ".szplt";

    PGrid *evenG;
    evenG = new PGrid(rank,size,0);

    pMat *loadMat, *evenMatFromLoad, *evenMat;


    tecIO *dataset1;
    vector<string> token;
    token.push_back("input");
    token.push_back("tecplot");
    token.push_back(prefix);
    token.push_back(suffix);
    token.push_back("150000");
    token.push_back("150100");
    token.push_back("10");
    token.push_back("Static_Pressure");
    token.push_back("-2");
    token.push_back("U");
    token.push_back("-1");
    token.push_back("V");
    token.push_back("-1");
   // token.push_back("CH4_mf");
   // token.push_back("1.0");
    dataset1 = new tecIO(token);
    string outdir="out";
    string outfile="stuff";


    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,1,0.0);
    dataset1->batchRead(loadMat);

    evenMat=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);

    evenMat->changeContext(loadMat);
        
    evenMatFromLoad = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    dataset1->batchRead(evenMatFromLoad);
    evenMatFromLoad->write_bin("tecplot.bin");

    assert((*evenMat)==(*evenMatFromLoad));   


    dataset1->calcAvg(evenMatFromLoad);
    dataset1->subAvg(evenMatFromLoad);

    dataset1->calcAvg(loadMat);
    dataset1->subAvg(loadMat);


    dataset1->calcNorm(loadMat);
    evenMat->changeContext(loadMat);
    assert((*evenMatFromLoad)==(*evenMatFromLoad));

    dataset1->normalize(loadMat);
    evenMat->changeContext(loadMat);


    for(int i=0;i<dataset1->numVars;i++)
    {
        dataset1->normFactor[i]=-1;
    }
    dataset1->normFactor[0]=-2;
    
    dataset1->calcNorm(evenMatFromLoad);
    dataset1->normalize(evenMatFromLoad);
    
    assert((*evenMat)==(*evenMatFromLoad));   


    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}