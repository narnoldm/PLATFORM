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
    if (rank != 0)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    cout << "testing tecplot metadata" << endl;

    

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
    token.push_back("151000");
    token.push_back("10");
    token.push_back("Static_Pressure");
    token.push_back("-1");
    token.push_back("Temperature");
    token.push_back("-2");
    //token.push_back("Temperature");
    //token.push_back("-1");
    dataset1 = new tecIO(token);
    string outdir="out2";
    string outfile="U";


    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,1,0.0);
    dataset1->batchRead(loadMat);


    dataset1->calcAvg(loadMat);
    dataset1->subAvg(loadMat);
    dataset1->calcNorm(loadMat);
    dataset1->normalize(loadMat);


    evenMat=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    evenMat->changeContext(loadMat);
    delete loadMat;

    pMat *U,*VT;
    vector<double> S;

    U=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    VT=new pMat(dataset1->nSets,dataset1->nSets,evenG,0,0,0.0);
    S.resize(dataset1->nSets);

    if(dataset1->nPoints/dataset1->nSets >=100)
    {
        evenMat->mos_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S);
    }
    else
    {
        evenMat->svd_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S);
    }

    printASCIIVecP0("S.txt",S.data(),S.size());
    VT->write_bin("VT.bin");


    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,1,0.0);
    loadMat->changeContext(U);

    tecIO *Uout=new tecIO();
    Uout->snap0 = 1;
    Uout->snapF = loadMat->N;
    Uout->snapSkip = 1;
    Uout->nSets = loadMat->N;
    Uout->prefix = "U";
    Uout->suffix = ".szplt";
    Uout->isInit = true;
    Uout->meshFile = dataset1->prefix+std::to_string(dataset1->snap0)+dataset1->suffix;
    Uout->fixedMesh = true;
    Uout->getDimNodes();
    Uout->varName=dataset1->varName;
    Uout->varIndex=dataset1->varIndex;
    Uout->numVars=Uout->varName.size();
    Uout->nPoints = Uout->nCells*Uout->numVars;


    Uout->activateGEMSbin("");
    Uout->batchWrite(loadMat);


    delete loadMat,evenMatFromLoad,evenMat;
    delete evenG,dataset1;


    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}