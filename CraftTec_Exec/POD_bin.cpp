#include "metadata.hpp"
#include "param.hpp"

using namespace ::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    paramMap inputFile("POD_bin.inp", rank);

    int debug_proc = 0;
    inputFile.getParamInt("stdout_proc", debug_proc);
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }


    int modeStart, modeEnd;
    inputFile.getParamInt("modeStart", modeStart);
    inputFile.getParamInt("modeEnd", modeEnd);

    string avgFile = "";
    bool readAvg = false;
    inputFile.getParamBool("readAvg", readAvg);
    if (readAvg)
        inputFile.getParamString("avgFile", avgFile);

    string input = "";
    inputFile.getParamString("inputString", input);
    cout << "input string is: " << input << endl;
    vector<string> token;
    tokenparse(input, "|", token);

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);

    pMat *evenMat;

    meta *dataset1;
    dataset1 = new tecIO(token);
    string outdir = "out2";
    string outfile = "U";

    evenMat = new pMat(dataset1->nPoints, dataset1->nSets, evenG, 0, 0, 0.0);
    dataset1->batchRead(evenMat);


    pMat *V, *YT;
    vector<double> S;

    int nModes = modeEnd - modeStart + 1;
    V = new pMat(dataset1->nPoints, nModes, evenG, 0, 0, 0.0);
    YT = new pMat(dataset1->nSets, dataset1->nSets, evenG, 0, 0, 0.0);

    S.resize(dataset1->nSets);

    evenMat->svd_run(dataset1->nPoints, dataset1->nSets, 0, 0, V, YT, S);

    printASCIIVecP0("S.txt", S, S.size());

    delete evenMat;

    YT->write_bin("YT.bin");

    meta *Uout = new tecIO();
    Uout->snap0 = 1;
    Uout->snapF = V->N;
    Uout->snapSkip = 1;
    Uout->nSets = V->N;
    Uout->prefix = "V";
    Uout->suffix = ".bin";
    Uout->isInit = true;
    Uout->nPoints = dataset1->nPoints;

    Uout->batchWrite(V, "Spatial_Modes", "Spatial_Mode_", modeStart - 1, modeEnd - 1, 1);

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}
