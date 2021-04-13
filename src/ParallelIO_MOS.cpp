#include "inputReader.hpp"

bool symprod(pMat *A, pMat *C, int n, int k)
{
    string UPLO = "U";
    string TRANS = "T";
    double ONE = 1.0;
    double ZERO = 0.0;
    int iONE = 1;

    pdsyrk(UPLO.c_str(), TRANS.c_str(), &n, &k, &ONE, A->dataD.data(), &iONE, &iONE, A->desc, &ZERO, C->dataD.data(), &iONE, &iONE, C->desc);
}

bool ParEig(pMat *corMat, double *eigVal, pMat *EigVec)
{
    int ibtype = 1;
    string jobz = "V";
    string range = "A";
    string UPLO = "U";
    string S = "S";
    double ONE = 1.0;
    double ZERO = 0.0;
    int iONE = 1;
    double orfac = 1e-3;
    double abstol = 2 * pdlamch(&(corMat->pG->icntxt), S.c_str());
    int m = 0;
    int nz = 0;

    vector<double> work;
    int lwork = -1;

    vector<int> iwork;
    int liwork = -1;

    vector<int> icluster;
    icluster.resize(2 * corMat->pG->prow * corMat->pG->pcol);

    vector<double> gap;
    gap.resize(corMat->pG->prow * corMat->pG->pcol);

    int info;

    vector<int> ifail;

    //pdsyevx(jobz.c_str(),range.c_str(),UPLO.c_str(),&(corMat->N),corMat->dataD.data(),&iONE,&iONE,corMat->desc,NULL,NULL,NULL,NULL,&abstol,&m,&nz,eigVal,&orfac,EigVec->dataD.data(),&iONE,&iONE,EigVec->desc,work.data(),&lwork,iwork.data(),&liwork,);
}

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //silence output for non root MPI processes
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();
    int debug_proc = 0;
    if (argc > 2)
    {
        debug_proc = atoi(argv[2]);
    }
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }
    string input = argv[1];
    cout << "input string is: " << input << endl;
    vector<string> tokens;
    tokenparse(input, "|", tokens);

    tecIO *dataset1 = new tecIO(tokens);
    //check for non standard input file name

    PGrid *loadGrid = new PGrid(rank, size, 1);
    pMat *data = new pMat(dataset1->nPoints, dataset1->nSets, loadGrid, 0, 1, 0.0);
    pMat *CorMat = new pMat(dataset1->nSets, dataset1->nSets, loadGrid, 0, 1, 0.0);

    dataset1->batchRead(data);
    dataset1->calcAvg(data);
    dataset1->subAvg(data);

    symprod(data, CorMat, dataset1->nSets, dataset1->nPoints);

    MPI_Barrier(MPI_COMM_WORLD);
    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}
