#include "metadata.hpp"
#include "param.hpp"

using namespace ::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    paramMap inputFile("POD_tec.inp", rank);

    int debug_proc = 0;
    inputFile.getParamInt("stdout_proc", debug_proc);
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    int mosStep = 0;
    inputFile.getParamInt("mosStep", mosStep);
    int modeStart, modeEnd;
    inputFile.getParamInt("modeStart", modeStart);
    inputFile.getParamInt("modeEnd", modeEnd);

    string avgFile = "";
    bool subAvg = true;
    bool readAvg = false;
    inputFile.getParamBool("subAvg", subAvg);
    if (subAvg) {
        inputFile.getParamBool("readAvg", readAvg);
        if (readAvg)
            inputFile.getParamString("avgFile", avgFile);
    }

    string input = "";
    inputFile.getParamString("inputString", input);
    cout << "input string is: " << input << endl;
    vector<string> token;
    tokenparse(input, "|", token);

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);

    pMat *evenMat;

    tecIO *dataset1;
    dataset1 = new tecIO(token);
    string outdir = "out2";
    string outfile = "U";

    if ((mosStep == 0) || (mosStep == 1) || (mosStep == 3))
    {
        evenMat = new pMat(dataset1->nPoints, dataset1->nSets, evenG, 0, 0, 0.0);
        dataset1->batchRead(evenMat);

        if (subAvg) {
            if (readAvg) {
                dataset1->readAvg(avgFile);
            } else {
                dataset1->calcAvg(evenMat);
            }
            dataset1->subAvg(evenMat);
        }
        dataset1->calcNorm(evenMat);
        dataset1->normalize(evenMat);
        //evenMat->write_bin("A.bin");
    }
    else
    {
    }

    pMat *U, *VT;
    vector<double> S;

    if ((mosStep == 0) || (mosStep == 3))
    {
        int nModes = modeEnd - modeStart + 1;
        U = new pMat(dataset1->nPoints, nModes, evenG, 0, 0, 0.0);
        VT = new pMat(dataset1->nSets, dataset1->nSets, evenG, 0, 0, 0.0);
    }
    if ((mosStep == 0) || (mosStep == 2) || (mosStep == 3))
    {
        S.resize(dataset1->nSets);
    }

    if ((dataset1->nPoints / dataset1->nSets >= 100) || (mosStep > 0))
    {
        if (mosStep == 0)
        {
            evenMat->mos_run(dataset1->nPoints, dataset1->nSets, 0, 0, U, VT, S);
        }
        else
        {
            evenMat->mos_run(dataset1->nPoints, dataset1->nSets, 0, 0, U, VT, S, modeStart, modeEnd, mosStep, evenG);
        }
    }
    else
    {
        evenMat->svd_run(dataset1->nPoints, dataset1->nSets, 0, 0, U, VT, S);
    }

    if ((mosStep == 0) || (mosStep == 2))
    {
        if (evenG->rank == 0)
        {
            printASCIIVecP0("S.txt", S, S.size());
        }
    }

    if ((mosStep == 0) || (mosStep == 1) || (mosStep == 3))
    {
        delete evenMat;
    }

    if ((mosStep == 0) || (mosStep == 3))
    {

        VT->write_bin("VT.bin");

        tecIO *Uout = new tecIO();
        Uout->snap0 = 1;
        Uout->snapF = U->N;
        Uout->snapSkip = 1;
        Uout->nSets = U->N;
        Uout->prefix = "U";
        Uout->suffix = ".szplt";
        Uout->isInit = true;
        Uout->meshFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;
        Uout->fixedMesh = true;
        Uout->getDimNodes();
        Uout->varName = dataset1->varName;
        Uout->varIndex = dataset1->varIndex;
        Uout->numVars = Uout->varName.size();
        Uout->nPoints = Uout->nCells * Uout->numVars;

        string firstFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;

        Uout->activateGEMSbin(firstFile.c_str());
        Uout->batchWrite(U, "Spatial_Modes", "Spatial_Mode_", modeStart - 1, modeEnd - 1, 1);
    }

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}
