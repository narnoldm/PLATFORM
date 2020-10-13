#include "metadata.hpp"

using namespace ::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
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
    vector<string> token;
    tokenparse(input, "|", token);

    tecIO *dataset1 = new tecIO(token);

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);

    pMat *A;

    A = new pMat(dataset1->nPoints, dataset1->nSets, evenG, 0, 0, 0.0);
    dataset1->batchRead(A);
    dataset1->calcAvg(A);
    dataset1->subAvg(A);
    dataset1->calcNorm(A);
    dataset1->normalize(A);

    int M=dataset1->nPoints, N=dataset1->nSets;


    int numModes = std::min(M,N);


    pMat *U, *VT, *UsT;
    vector<double> S;

    U = new pMat(M, std::min(M, N), evenG);
    VT = new pMat(std::min(M, N), N, evenG);
    S.resize(std::min(M, N));

    A->svd_run(M, N, 0, 0, U, VT, S);
    for (int i = 0; i < numModes; i++)
        cout << S[i] << endl;
    vector<int> P;

    UsT = new pMat(numModes, U->M, evenG);
    UsT->transpose(U, UsT->M, UsT->N, 0, 0);

    UsT->qr_run(UsT->M, UsT->N, 0, 0, P);

    vector<int> gP;
    if (!rank)
    {
        readMat("P.bin", gP);
        gP.resize(numModes);
        for (int i = 0; i < gP.size(); i++)
        {
            gP[i]--; //switch to 0 indexing
            //switch to physical points
            gP[i]=gP[i]%dataset1->nCells;



            
            cout << gP[i] << endl;
        }
    }











    string filename = "gemsma1.bin";

    if (A->check_bin_size(filename, M, N))
    {

        A = new pMat(M, N, evenG);
        A->read_bin(filename);

        pMat *U, *VT, *UsT;
        vector<double> S;

        U = new pMat(M, std::min(M, N), evenG);
        VT = new pMat(std::min(M, N), N, evenG);
        S.resize(std::min(M, N));

        A->svd_run(M, N, 0, 0, U, VT, S);
        for (int i = 0; i < numModes; i++)
            cout << S[i] << endl;
        vector<int> P;

        UsT = new pMat(numModes, U->M, evenG);
        UsT->transpose(U, UsT->M, UsT->N, 0, 0);

        UsT->qr_run(UsT->M, UsT->N, 0, 0, P);

        vector<int> gP;
        if (!rank)
        {
            readMat("P.bin", gP);
            gP.resize(numModes);
            for (int i = 0; i < gP.size(); i++)
            {
                gP[i]--; //switch to 0 indexing
                cout << gP[i] << endl;
                //switch to physical points
            }
            //add boundary points

            //remove duplicate physical points
        }
    }

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}
