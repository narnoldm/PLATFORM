

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

    cout << "testing binaryset metadata" << endl;

    

    string prefix = "testsplit/test";
    string suffix = ".bin";

    PGrid *evenG;
    evenG = new PGrid(rank,size,0);

    pMat *loadMat, *evenMatFromLoad, *evenMat;

    string outdir="out";
    string outfile="U";

    meta *dataset1;
    vector<string> token;
    token.push_back("input");
    token.push_back("binaryset");
    token.push_back(prefix);
    token.push_back(suffix);
    token.push_back("1");
    token.push_back("10");
    token.push_back("1");
    //token.push_back("Temperature");
    //token.push_back("-1");
    dataset1 = new meta(token);


    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,1,0.0);
    dataset1->batchRead(loadMat);
    evenMat=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    evenMat->changeContext(loadMat);

    pMat *U,*VT;
    vector<double> S;

    U=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    VT=new pMat(dataset1->nSets,dataset1->nSets,evenG,0,0,0.0);
    S.resize(dataset1->nSets);
    evenMat->svd_run(dataset1->nPoints,dataset1->nSets,0,0,U,VT,S);
    evenMat->outerProductSum(U,'N',VT,'N',S,0);

    for(int i=0;i<S.size();i++)
        cout<<S[i]<<endl;


    pMat *U2,*VT2;
    vector<double> S2;

    U2=new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,0,0.0);
    VT2=new pMat(dataset1->nSets,dataset1->nSets,evenG,0,0,0.0);
    S2.resize(dataset1->nSets);

    evenMat->mos_run(dataset1->nPoints,dataset1->nSets,0,0,U2,VT2,S2);

    for(int i=0;i<S2.size();i++)
        cout<<S2[i]<<endl;

    loadMat = new pMat(dataset1->nPoints,dataset1->nSets,evenG,0,1,0.0);
    loadMat->changeContext(U);
    dataset1->batchWrite(loadMat,outdir,outfile);

    delete loadMat,evenMatFromLoad,evenMat;
    delete evenG,dataset1;


    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}