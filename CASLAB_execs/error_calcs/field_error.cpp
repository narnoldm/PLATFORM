#include "metadata.hpp"
#include "param.hpp"

using namespace :: std;

int main(int argc, char *argv[]) {

    // some setup
	MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    // must provide input file as command-line argument
    assert (argc == 2);
    string inpFile = argv[1];
    paramMap inputFile(inpFile);

    int debug_proc;
    inputFile.getParamInt("debug_proc", debug_proc, 0);
    if (rank != debug_proc) {
        std::cout.rdbuf(sink.rdbuf());
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
        tecIO* setROM = new tecIO(token);

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
        errDir = romDir + "/k" + to_string(setBasis->nSets);
        errSuffix = to_string(setROM->snap0) + "_" + to_string(setROM->snapF) + "_" + to_string(setROM->snapSkip);
    } else {
        errDir = fomDir + "/k" + to_string(setBasis->nSets);
        errSuffix = to_string(setFOM->snap0) + "_" + to_string(setFOM->snapF) + "_" + to_string(setFOM->snapSkip);
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
        setBasis->batchRead(basis);

        // project data to low-dimension
        latentCode = new pMat(setBasis->nSets, setFOM->nSets, evenG, false);
        QTruth_proj = new pMat(setFOM->nPoints, setFOM->nSets, evenG, false);
        latentCode->matrix_Product('T', 'N', basis->N, QTruth->N, basis->M, basis, 0, 0, QTruth, 0, 0, 1.0, 0.0, 0, 0);
        QTruth_proj->matrix_Product('N', 'N', basis->M, QTruth->M, basis->N, basis, 0, 0, latentCode, 0, 0, 1.0, 0.0, 0, 0);

        // de-normalize and de-center, if normalization/centering was requested
        if (normalize)
            setFOM->unNormalize(QTruth_proj);
        if (center)
            setFOM->addAvg(QTruth_proj);

        if (errType == 1) {
            QComp = QTruth_proj;
        } else {
            // don't need original FOM solution any more
            delete QTruth;
            QTruth = QTruth_proj;
        }

        // output projected field or latent code history, if requested
        if (outProjField)
            setFOM->batchWrite(QTruth_proj, projDir, "proj_sol_");
        if (outLatentCode)
            latentCode->write_bin(projDir + "/latentCode.bin");

        delete latentCode;

    } else {
        QComp = new pMat(setROM->nPoints, setROM->nSets, evenG, false);
    }

    // load ROM data
    if (errType > 1) {
        setROM->batchRead(QComp);
    }

    // ----- FINISH DATA LOADING AND PREPROCESSING -----

    // ----- COMPUTE ERROR, WRITE OUTPUTS -----

    pMat* ones = new pMat(1, setFOM->nCells, evenG, 0, 0, 1.0, false);
    pMat* errField = new pMat(setFOM->nPoints, setFOM->nSets, evenG, false);
    pMat* errVar = new pMat(setFOM->numVars, setFOM->nSets, 0, 0, 0.0, false);
    pMat* norm = new pMat(setFOM->numVars, setFOM->nSets, 0, 0, 0.0, false);

    // # absolute error #
    // full field error snapshots [abs(qTruth - qComp)]
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(QTruth->dataD[i] - QComp->dataD[i]);
    if (outAbsErrField)
        setFOM->batchWrite(errField, errDir, "abs_err_");

    // variable error history [sum(abs(qTruth - qComp)) / nCells]
    for (int i = 0; i < setFOM->numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, ones, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= setFOM->nCells;
    errVar->write_bin(errDir + "/abs_avg_err_" + errSuffix + ".bin");

    // relative variable error history [sum(abs(qTruth - qComp)) / sum(abs(qTruth))]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(QTruth->dataD[i]);
    for (int i = 0; i < setFOM->numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, ones, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= norm->dataD[i];
    errVar->write_bin(errDir + "/abs_rel_err_" + errSuffix + ".bin");

    // # L2 error #
    for (int i = 0; i < errField->dataD.size(); ++i) {
        errField->dataD[i] = QTruth->dataD[i] - QComp->dataD[i];
        errField->dataD[i] = errField->dataD[i] * errField->dataD[i];
    }

    // variable error history [||qTruth - qComp|| / nCells]
    for (int i = 0; i < setFOM->numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, ones, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] = sqrt(errVar->dataD[i]) / setFOM->nCells;
    errVar->write_bin(errDir + "/l2_avg_err_" + errSuffix + ".bin");

    // relative variable error history [||qTruth - qComp|| / ||qTruth||]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = QTruth->dataD[i] * QTruth->dataD[i];
    for (int i = 0; i < setFOM->numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setFOM->nSets, setFOM->nCells, ones, 0, 0, errField, i * setFOM->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= sqrt(norm->dataD[i]);
    errVar->write_bin(errDir + "/l2_rel_err_" + errSuffix + ".bin");

    // ----- FINISH COMPUTE ERROR, WRITE OUTPUTS -----

    // cleanup
    delete QTruth_proj;
    delete QComp;
    delete ones;
    delete errField;
    delete errVar;

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}