#include "metadata.hpp"
#include "param.hpp"
#include "extern_func.hpp"

using namespace :: std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    paramMap inputFile("POD_tec.inp", rank);

    int debug_proc = 0;
    inputFile.getParamInt("stdout_proc", debug_proc);
    if (rank != debug_proc)
    {
        cout.rdbuf(sink.rdbuf());
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

    string cellIDFile;
    inputFile.getParamString("cellIDFile", cellIDFile, "");

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);

    pMat *evenMat;

    tecIO *dataset1;
    dataset1 = new tecIO(token, cellIDFile);

    double t1;

    if ((mosStep == 0) || (mosStep == 1) || (mosStep == 3))
    {

        evenMat = new pMat(dataset1->nPoints, dataset1->nSets, evenG, 0, 0, 0.0, false);
        dataset1->activateReorder();
        dataset1->batchRead(evenMat);

        // read centering inputs
        string centerFile, centerMethod;
        bool center, centerIsField, writeCentering;
        inputFile.getParamBool("center", center);
        if (center)
        {
            inputFile.getParamString("centerFile", centerFile, "");
            inputFile.getParamString("centerMethod", centerMethod, "");
            inputFile.getParamBool("centerIsField", centerIsField, false);
            inputFile.getParamBool("writeCentering", writeCentering, true);
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

            t1 = MPI_Wtime();
            if (centerFile != "")
            {
                dataset1->calcCentering(evenMat, centerFile, true, writeCentering);
            }
            else
            {
                dataset1->calcCentering(evenMat, centerMethod, centerIsField, writeCentering);
            }
            dataset1->centerData(evenMat, false);
            cout << "Centering time: " << MPI_Wtime() - t1 << endl;
        }

        // read scaling inputs
        string scaleFile, scaleMethod;
        bool scale, scaleIsField, writeScaling;
        inputFile.getParamBool("scale", scale);
        if (scale)
        {
            inputFile.getParamString("scaleFile", scaleFile, "");
            inputFile.getParamString("scaleMethod", scaleMethod, "");
            inputFile.getParamBool("scaleIsField", scaleIsField, false);
            inputFile.getParamBool("writeScaling", writeScaling, true);
            if ((scaleFile == "") && (scaleMethod == ""))
            {
                cout << "Must provide scaleFile or scaleMethod if scale = true" << endl;
                throw(-1);
            }
            if ((scaleFile != "") && (scaleMethod != ""))
            {
                cout << "Can only set scaleFile OR scaleMethod if center = true" << endl;
                throw(-1);
            }

            t1 = MPI_Wtime();
            if (scaleFile != "")
            {
                dataset1->calcScaling(evenMat, scaleFile, true, writeScaling);
            }
            else
            {
                dataset1->calcScaling(evenMat, scaleMethod, scaleIsField, writeScaling);
            }
            dataset1->scaleData(evenMat, false);
            cout << "Scaling time: " << MPI_Wtime() - t1 << endl;
        }

    }

    pMat *U, *VT;
    vector<double> S;

    if ((mosStep == 0) || (mosStep == 3))
    {
        int nModes = modeEnd - modeStart + 1;
        U = new pMat(dataset1->nPoints, nModes, evenG, 0, 0, 0.0, false);
        VT = new pMat(dataset1->nSets, dataset1->nSets, evenG, 0, 0, 0.0, false);
    }
    if ((mosStep == 0) || (mosStep == 2) || (mosStep == 3))
    {
        S.resize(dataset1->nSets);
    }

    if (mosStep == 0)
    {
        evenMat->mos_run(dataset1->nPoints, dataset1->nSets, 0, 0, U, VT, S, modeStart, modeEnd);
    }
    else
    {
        evenMat->mos_run(dataset1->nPoints, dataset1->nSets, 0, 0, U, VT, S, modeStart, modeEnd, mosStep, evenG);
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

        // determine output format
        int outMode;
        bool writeModesSZPLT, writeModesBin;
        inputFile.getParamBool("writeModesSZPLT", writeModesSZPLT, true);
        inputFile.getParamBool("writeModesBin", writeModesBin, true);

        // set up meta
        tecIO *Uout = new tecIO();
        Uout->snap0 = 1;
        Uout->snapF = U->N;
        Uout->snapSkip = 1;
        Uout->nSets = U->N;
        Uout->prefix = "U";
        Uout->suffix = ".szplt";
        Uout->isInit = true;
        Uout->meshFile = dataset1->prefix + to_string(dataset1->snap0) + dataset1->suffix;
        Uout->fixedMesh = true;
        Uout->getDimNodes();
        Uout->varName = dataset1->varName;
        Uout->varIndex = dataset1->varIndex;
        Uout->numVars = Uout->varName.size();
        Uout->nPoints = Uout->nCells * Uout->numVars;
        Uout->cellIDFile = dataset1->cellIDFile;

        // for reordering output
        Uout->activateReorder();

        if (writeModesSZPLT)
        {
            if (writeModesBin)
            {
                outMode = 0;
            }
            else
            {
                outMode = 1;
            }

            // write modes to file
            Uout->batchWrite(U, "Spatial_Modes", "Spatial_Mode_", 0, modeEnd - modeStart + 1, 1, modeStart, 1, 0, outMode);

        }
        else if (writeModesBin)
        {
            Uout->batchWrite_bin(U, "Spatial_Modes", "Spatial_Mode_", 0, modeEnd - modeStart + 1, 1, modeStart, 1);
        }
        else
        {
            cout << "No output requested" << endl;
        }

    }

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}
