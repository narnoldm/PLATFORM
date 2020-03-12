#include "inputReader.hpp"


int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    //silence output for non root MPI processes
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
    vector<string> tokens;
    tokenparse(input,"|",tokens);

    tecIO *dataset1=new tecIO(tokens);
    //check for non standard input file name

    PGrid *loadGrid = new PGrid(rank,size,1);
    pMat *data = new pMat(dataset1->nPoints,dataset1->nSets,loadGrid,0,1,0.0);

    dataset1->batchRead(data);
    dataset1->calcAvg(data);
    dataset1->subAvg(data);
    dataset1->calcNorm(data);
    dataset1->normalize(data);

    PGrid *evenGrid = new PGrid(rank,size,0);
    pMat *dataE = new pMat(dataset1->nPoints,dataset1->nSets,evenGrid,0,0,0.0);
    dataE->changeContext(data);
    delete loadGrid;
    delete data;

    pMat *U=new pMat(dataset1->nPoints,dataset1->nSets,evenGrid,0,0,0.0);
    pMat *VT=new pMat(dataset1->nSets,dataset1->nSets,evenGrid,0,0,0.0);
    
    vector<double> S;
    S.resize(dataset1->nSets);

    dataE->svd_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S);


    VT->write_bin("VT.bin");
    printASCIIVecP0("S.txt",S.data(),S.size());
    delete VT;
    S.clear();

    loadGrid = new PGrid(rank,size,1);
    data = new pMat(dataset1->nPoints,dataset1->nSets,loadGrid,0,1,0.0);

    data->changeContext(dataE);
    delete dataE,evenGrid;

    dataset1->batchWrite(data,"Spatial_Modes","Spatial_Mode");
    

    
    MPI_Barrier(MPI_COMM_WORLD);
    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}





