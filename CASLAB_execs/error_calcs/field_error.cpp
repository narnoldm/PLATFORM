#include "metadata.hpp"
#include "param.hpp"

using namespace :: std;

void calc_integrated_error(pMat* dataMat, int rank, vector<string>& varNames, string outFile) {
    // errDir + "/abs_avg_int_err" + errSuffix + ".dat"

    pMat* onesCol = new pMat(dataMat->N, 1, dataMat->pG, 0, 0, 1.0, false);
    pMat* dataMatInt = new pMat(dataMat->M, 1, dataMat->pG, false);
    pMat* dataMatIntP0 = new pMat(dataMat->M, 1, dataMat->pG, 0, 2, 0.0, false);
    double dataMatIntAvg = 0.0;

    dataMatInt->matrix_Product('N', 'N', dataMat->M, 1, dataMat->N, dataMat, 0, 0, onesCol, 0, 0, 1.0, 0.0, 0, 0);
    dataMatIntP0->changeContext(dataMatInt, false);
    if (!rank) {
        ofstream out;
		out.open(outFile, ios::trunc);
		for (int i = 0; i < dataMatIntP0->dataD.size(); ++i) {
            out << varNames[i] + ": " << setprecision(numeric_limits<double>::digits10) << dataMatIntP0->dataD[i] << endl;
            dataMatIntAvg += dataMatIntP0->dataD[i];
        }
        dataMatIntAvg /= dataMat->M;
        out << "Average: " << setprecision(numeric_limits<double>::digits10) << dataMatIntAvg << endl;
        out.close();
    }
}

int main(int argc, char *argv[]) {

    // some setup
	MPI_Init(&argc, &argv);
    int rank, size, ierr;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    // must provide input file as command-line argument
    assert (argc == 2);
    string inpFile = argv[1];
    paramMap inputFile(inpFile);

    int debug_proc = 0;
    inputFile.getParamInt("debug_proc", debug_proc, 0);
    if (rank != debug_proc) {
        cout.rdbuf(sink.rdbuf());
    }

    // ----- INPUT PARAMETERS -----

    string fomInputString, romInputString, basisInputString, centerFile, normFile;
    bool center, normalize, outProjField, outLatentCode, outAbsErrField;
    int errType;

    // 1: compute error between FOM and projected FOM solution
    // 2: compute error between FOM and ROM solution
    // 3: compute error between projected FOM and ROM solution
    inputFile.getParamInt("errType", errType);

    inputFile.getParamString("fomInputString", fomInputString);  // input token for FOM data series

    // input token for ROM data series
    if (errType > 1)
        inputFile.getParamString("romInputString", romInputString);

    // if FOM projection is required
    bool projFOM = false;
    if ((errType == 1) || (errType == 3))
        projFOM = true;

    if (projFOM) {
        // input token for spatial mode series
        inputFile.getParamString("basisInputString", basisInputString);

        // centering FOM data before projection
        inputFile.getParamBool("center", center);
        if (center) {
            // path to data centering profile
            // if not provided, use mean field
            inputFile.getParamString("centerFile", centerFile, "");
        }

        // normalizing FOM data before projection (after centering, if requested)
        inputFile.getParamBool("normalize", normalize);
        if (normalize) {
            // path to data normalization profile
            // if not provided, use normalization constants provided in fomInputString
            inputFile.getParamString("normFile", normFile, "");
        }

        inputFile.getParamBool("outProjField", outProjField, false);  // output unsteady projected FOM solutions
        inputFile.getParamBool("outLatentCode", outLatentCode, false);  // output unsteady low-dimensional state

    }

    inputFile.getParamBool("outAbsErrField", outAbsErrField, false);  // output absolute unsteady error fields

    // ----- FINISH INPUT PARAMETERS -----

    // ----- PREPROCESSING -----

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);
    vector<string> token;
    string fomDir, romDir, projDir, errDir, errSuffix;

    // FOM dataset
    tokenparse(fomInputString, "|", token);
    tecIO* setFOM = new tecIO(token);

    // ROM dataset
    tecIO* setROM;
    if (errType > 1) {
        token.clear();
        tokenparse(romInputString, "|", token);
        setROM = new tecIO(token);

        // error checking
        assert (setFOM->dim == setROM->dim);
        assert (setFOM->nCells == setROM->nCells);
        assert (setFOM->numVars == setROM->numVars);
        assert (setFOM->nSets == setROM->nSets);

    }

    // basis dataset
    meta* setBasis;
    if (projFOM) {
        token.clear();
        tokenparse(basisInputString, "|", token);
        setBasis = new meta(token);

        // error checking
        assert (setBasis->nPoints == setFOM->nPoints);

    }

    // output directories
    size_t fomDirPos = setFOM->prefix.find_last_of("/");
    fomDir = setFOM->prefix.substr(0, fomDirPos);
    if (projFOM) {
        projDir = fomDir + "/k" + to_string(setBasis->nSets);
        if (outProjField || outLatentCode) {
            if (!rank)
                system(("mkdir " + projDir).c_str());
        }
    }

    if (errType > 1) {
        size_t romDirPos = setROM->prefix.find_last_of("/");
        romDir = setROM->prefix.substr(0, romDirPos);
        errDir = romDir;
        errSuffix = "_" + to_string(setROM->snap0) + "_" + to_string(setROM->snapF) + "_" + to_string(setROM->snapSkip);
        if (errType == 2) {
            errSuffix = "_vs_raw" + errSuffix;
        } else {
            errSuffix = "_vs_proj" + errSuffix;
            errDir += "/k" + to_string(setBasis->nSets);
        }
    } else {
        errDir = fomDir + "/k" + to_string(setBasis->nSets);
        errSuffix = "_" + to_string(setFOM->snap0) + "_" + to_string(setFOM->snapF) + "_" + to_string(setFOM->snapSkip);
    }
    if (!rank)
        system(("mkdir " + errDir).c_str());

    // ----- FINISH PREPROCESSING -----

    // ----- DATA LOADING -----

    // load FOM data
    pMat* QTruth = new pMat(setFOM->nPoints, setFOM->nSets, evenG, false);
    setFOM->batchRead(QTruth);

    pMat *QComp, *QTruth_proj, *latentCode, *basis;
    // project FOM data if requested
    if (projFOM) {

        // center data, if requested
        if (center) {
            if (centerFile == "") {
                setFOM->calcAvg(QTruth);
            } else {
                setFOM->readAvg(centerFile);
            }
            setFOM->subAvg(QTruth);
        }

        // normalized data, if requested
        if (normalize) {
            if (normFile == "") {
                setFOM->calcNorm(QTruth);
            } else {
                cout << "Norm file read not implemented yet" << endl;
                throw(-1);
            }
            setFOM->normalize(QTruth);
        }

        // load trial basis
        pMat* basis = new pMat(setBasis->nPoints, setBasis->nSets, evenG, false);
        setBasis->batchRead(basis);

        // project data to low-dimension
        latentCode = new pMat(setBasis->nSets, setFOM->nSets, evenG, false);
        QTruth_proj = new pMat(setFOM->nPoints, setFOM->nSets, evenG, false);
        latentCode->matrix_Product('T', 'N', basis->N, QTruth->N, basis->M, basis, 0, 0, QTruth, 0, 0, 1.0, 0.0, 0, 0);
        QTruth_proj->matrix_Product('N', 'N', basis->M, latentCode->N, basis->N, basis, 0, 0, latentCode, 0, 0, 1.0, 0.0, 0, 0);

        // de-normalize and de-center, if normalization/centering was requested
        if (normalize) {
            setFOM->unNormalize(QTruth_proj);
            if (errType == 1)
                setFOM->unNormalize(QTruth);
        }
        if (center) {
            setFOM->addAvg(QTruth_proj);
            if (errType == 1)
                setFOM->addAvg(QTruth);
        }

        if (errType == 1) {
            QComp = QTruth_proj;
        } else {
            // don't need original FOM solution any more
            destroyPMat(QTruth, false);
            QTruth = QTruth_proj;
        }

        // output projected field or latent code history, if requested
        if (outProjField)
            setFOM->batchWrite(QTruth_proj, projDir, "proj_sol_");
        if (outLatentCode)
            latentCode->write_bin(projDir + "/latentCode.bin");

        destroyPMat(latentCode, false);

    }

    if (errType > 1)
        QComp = new pMat(setROM->nPoints, setROM->nSets, evenG, false);

    // load ROM data
    if (errType > 1) {
        setROM->batchRead(QComp);
    }

    // ----- FINISH DATA LOADING AND PREPROCESSING -----

    // ----- COMPUTE ERROR, WRITE OUTPUTS -----

    pMat* onesRow = new pMat(1, setFOM->nCells, evenG, 0, 0, 1.0, false);
    pMat* errField = new pMat(setFOM->nPoints, setFOM->nSets, evenG, false);
    pMat* errVar = new pMat(setFOM->numVars, setFOM->nSets, evenG, false);
    pMat* norm = new pMat(setFOM->numVars, setFOM->nSets, evenG, false);

    // # absolute error #
    // full field error snapshots [abs(qTruth - qComp)]
    cout << "Calculating absolute error" << endl;
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(QTruth->dataD[i] - QComp->dataD[i]);
    if (outAbsErrField)
        setFOM->batchWrite(errField, errDir, "abs_err_");

    // variable error history [sum(abs(qTruth - qComp)) / nCells]
    for (int i = 0; i < setFOM->numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, onesRow, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= setFOM->nCells;
    errVar->write_bin(errDir + "/abs_avg_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, rank, setFOM->varName, errDir + "/abs_avg_sum_err" + errSuffix + ".dat");

    // relative variable error history [sum(abs(qTruth - qComp)) / sum(abs(qTruth))]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(QTruth->dataD[i]);
    for (int i = 0; i < setFOM->numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, onesRow, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= norm->dataD[i];
    errVar->write_bin(errDir + "/abs_rel_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, rank, setFOM->varName, errDir + "/abs_rel_sum_err" + errSuffix + ".dat");

    // # L2 error #
    cout << "Calculating L2 error" << endl;
    for (int i = 0; i < errField->dataD.size(); ++i) {
        errField->dataD[i] = QTruth->dataD[i] - QComp->dataD[i];
        errField->dataD[i] = errField->dataD[i] * errField->dataD[i];
    }

    // variable error history [||qTruth - qComp|| / nCells]
    for (int i = 0; i < setFOM->numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, onesRow, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] = sqrt(errVar->dataD[i]) / setFOM->nCells;
    errVar->write_bin(errDir + "/l2_avg_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, rank, setFOM->varName, errDir + "/l2_avg_sum_err" + errSuffix + ".dat");

    // relative variable error history [||qTruth - qComp|| / ||qTruth||]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = QTruth->dataD[i] * QTruth->dataD[i];
    for (int i = 0; i < setFOM->numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, onesRow, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= sqrt(norm->dataD[i]);
    errVar->write_bin(errDir + "/l2_rel_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, rank, setFOM->varName, errDir + "/l2_rel_sum_err" + errSuffix + ".dat");

    // ----- FINISH COMPUTE ERROR, WRITE OUTPUTS -----

    // cleanup
    MPI_Barrier(MPI_COMM_WORLD);
    destroyPMat(QTruth, false);
    destroyPMat(QComp, false);
    destroyPMat(onesRow, false);
    destroyPMat(errField, false);
    destroyPMat(errVar, false);
    destroyPMat(norm, false);

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}