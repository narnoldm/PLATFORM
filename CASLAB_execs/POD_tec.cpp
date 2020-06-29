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
    int mosStep = 0;
    int modeStart, modeEnd;
    string avgFile;
    if(argc>2)
    {
        debug_proc=atoi(argv[2]);
        if(argc>3)
        {
            avgFile=argv[3];
             
            if (argc > 4) {
                mosStep = atoi(argv[4]);

                // modeStart and modeEnd
                if (argc > 5) {
                    modeStart = atoi(argv[5]);
                    if (argc > 6) {
                        modeEnd = atoi(argv[6]);
                    }
                }

            }

        }
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

    if ( (mosStep == 0) || (mosStep == 1) || (mosStep == 3) ) {
        evenMat=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
        dataset1->batchRead(evenMat);

        if(argc>3)
        {
            dataset1->readAvg(avgFile);
        }
        else
        {
            dataset1->calcAvg(evenMat);
        }
        dataset1->subAvg(evenMat);
        dataset1->calcNorm(evenMat);
        dataset1->normalize(evenMat);
        evenMat->write_bin("A.bin");
    } else {

    }
    
    pMat *U,*VT;
    vector<double> S;

    if ( (mosStep == 0) || (mosStep == 3) ) {
        if (argc < 7) {
            modeEnd = dataset1->nSets;
            if (argc < 6) {
                modeStart = 1;
            }
        }
        int nModes = modeEnd - modeStart + 1;
        U=new pMat(dataset1->nPoints,nModes,evenG,0,0,0.0);
        VT=new pMat(dataset1->nSets,dataset1->nSets,evenG,0,0,0.0);
    }
    if ( (mosStep == 0) || (mosStep == 2) || (mosStep == 3) ) {
        S.resize(dataset1->nSets);
    }

    if(dataset1->nPoints/dataset1->nSets >=100)
    {
        if (mosStep == 0) {
            evenMat->mos_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S);
        } else {
            evenMat->mos_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S,modeStart,modeEnd,mosStep,evenG);
        }
        
    }
    else
    {
        evenMat->svd_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S);
    }

    if ( (mosStep == 0) || (mosStep == 2) ) {
        if (evenG->rank == 0) {
            printASCIIVecP0("S.txt",S.data(),S.size());
        }
    }
    
    if ( (mosStep == 0) || (mosStep == 1) || (mosStep == 3) ) {
        delete evenMat;
    }
    
    if ( (mosStep == 0) || (mosStep == 3) ) {

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

        string firstFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;
        U->write_bin("U.bin");
        Uout->activateGEMSbin(firstFile.c_str());
        Uout->batchWrite(U,"Spatial_Modes", "Spatial_Mode_", modeStart-1, modeEnd-1, 1);

    }
    
    




    return 0;
}
