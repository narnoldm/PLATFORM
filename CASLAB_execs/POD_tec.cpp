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

    if ((mosStep == 0) || (mosStep == 1) || (mosStep == 3))
    {

        evenMat = new pMat(dataset1->nPoints, dataset1->nSets, evenG, 0, 0, 0.0);
        dataset1->batchRead(evenMat);

        // read centering inputs
        string centerFile, centerMethod;
        bool center, centerIsField;
        inputFile.getParamBool("center", center);
        if (center)
        {
            inputFile.getParamString("centerFile", centerFile, "");
            inputFile.getParamString("centerMethod", centerMethod, "");
            inputFile.getParamBool("centerIsField", centerIsField, false);
            if ((centerFile == "") && (centerMethod == ""))
            {
                cout << "Must provide centerFile or centerMethod if center = true" << endl;
                throw(-1);
            }
            if ((centerFile != "") && (centerMethod != ""))
            {
                cout << "Can only set centerFile OR centerMethod if center = true" << endl;
                throw(-1);
            }
                
            if (centerFile != "")
            {
                dataset1->calcCentering(evenMat, centerFile, true);
            }
            else
            {
                dataset1->calcCentering(evenMat, centerMethod, centerIsField);
            }
            dataset1->centerData(evenMat, false);
        }

        // read scaling inputs
        string scaleFile, scaleMethod;
        bool scale, scaleIsField;
        inputFile.getParamBool("scale", scale);
        if (scale)
        {
            inputFile.getParamString("scaleFile", scaleFile, "");
            inputFile.getParamString("scaleMethod", scaleMethod, "");
            inputFile.getParamBool("scaleIsField", scaleIsField, false);
            if ((scaleFile == "") && (scaleMethod == ""))
            {
                cout << "Must provide centerFile or centerMethod if center = true" << endl;
                throw(-1);
            }
            if ((scaleFile != "") && (scaleMethod != ""))
            {
                cout << "Can only set centerFile OR centerMethod if center = true" << endl;
                throw(-1);
            }
            if (scaleFile != "")
            {
                dataset1->calcScaling(evenMat, scaleFile, true);
            }
            else
            {
                dataset1->calcScaling(evenMat, scaleMethod, scaleIsField);
            }
            dataset1->scaleData(evenMat, false);
        }

        MPI_Barrier(MPI_COMM_WORLD);
        int ierr;
        MPI_Abort(MPI_COMM_WORLD, ierr);

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
        Uout->batchWrite(U, "Spatial_Modes", "Spatial_Mode_", 0, modeEnd - modeStart + 1, 1, modeStart, 1);
    }

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}
