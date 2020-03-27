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
    int debug_proc=0;
    if(argc>2)
    {
        debug_proc=atoi(argv[2]);
    }
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }
    string input = argv[1];
    cout<<"input string is: "<<input<<endl;
    vector<string> token;
    tokenparse(input,"|",token);



    PGrid *evenG;
    evenG = new PGrid(rank,size,0);

    pMat *evenMat;


    tecIO *dataset1;
    dataset1 = new tecIO(token);
    string outdir="out2";
    string outfile="U";


    evenMat=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    dataset1->batchRead(evenMat);

    
    dataset1->calcAvg(evenMat);
    dataset1->subAvg(evenMat);
    dataset1->calcNorm(evenMat);
    dataset1->normalize(evenMat);


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
    

    tecIO *Uout=new tecIO();
    Uout->snap0 = 1;
    Uout->snapF = U->N;
    Uout->snapSkip = 1;
    Uout->nSets = U->N;
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
    Uout->batchWrite(U);



    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}