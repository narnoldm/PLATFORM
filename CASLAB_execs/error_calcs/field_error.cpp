#include "metadata.hpp"
#include "param.hpp"
#include "error_funcs.hpp"

using namespace :: std;

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
    inputFile.getParamInt("stdout_proc", debug_proc, 0);
    if (rank != debug_proc) {
        cout.rdbuf(sink.rdbuf());
    }

    // ----- INPUT PARAMETERS -----

    string fomInputString, romInputString, basisInputString;
    string centerFile, centerMethod, scaleFile, scaleMethod;
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
            inputFile.getParamString("centerMethod", centerMethod, "");
        }

        // normalizing FOM data before projection (after centering, if requested)
        inputFile.getParamBool("normalize", normalize);
        if (normalize) {
            // path to data normalization profile
            // if not provided, use normalization constants provided in fomInputString
            inputFile.getParamString("scaleFile", scaleFile, "");
            inputFile.getParamString("scaleMethod", scaleMethod, "");
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
    setFOM->activateReorder(setFOM->prefix + to_string(setFOM->snap0) + setFOM->suffix);

    // ROM dataset
    tecIO* setROM;
    if (errType > 1) {
        token.clear();
        tokenparse(romInputString, "|", token);
        setROM = new tecIO(token);
        setROM->activateReorder(setROM->prefix + to_string(setROM->snap0) + setROM->suffix);

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
        projDir = fomDir + "/projection/k" + to_string(setBasis->nSets);
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
            errDir += "/projection/k" + to_string(setBasis->nSets);
        }
    } else {
        errDir = fomDir + "/projection/k" + to_string(setBasis->nSets);
        errSuffix = "_" + to_string(setFOM->snap0) + "_" + to_string(setFOM->snapF) + "_" + to_string(setFOM->snapSkip);
    }
    if (!rank)
        system(("mkdir -p " + errDir).c_str());

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
                setFOM->calcCentering(QTruth, centerMethod);
            } else {
                setFOM->readCentering(centerFile);
            }
            setFOM->centerData(QTruth);
        }

        // normalized data, if requested
        if (normalize) {
            if (scaleFile == "") {
                setFOM->calcScaling(QTruth, scaleMethod);
            } else {
                cout << "Norm file read not implemented yet" << endl;
                throw(-1);
            }
            setFOM->scaleData(QTruth);
        }
        
        if (outProjField)
            setFOM->batchWrite(QTruth, projDir, "fom_sol_raw_");

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
            setFOM->scaleData(QTruth_proj, true);
            if (errType == 1)
                setFOM->scaleData(QTruth, true);
        }
        if (center) {
            setFOM->centerData(QTruth_proj, true);
            if (errType == 1)
                setFOM->centerData(QTruth, true);
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

    // ----- FINISH DATA LOADING -----

    // ----- COMPUTE ERROR, WRITE OUTPUTS -----

    calc_abs_and_l2_error(QTruth, QComp, setFOM, errDir, errSuffix, outAbsErrField);

    // ----- FINISH COMPUTE ERROR, WRITE OUTPUTS -----

    // cleanup
    MPI_Barrier(MPI_COMM_WORLD);
    destroyPMat(QTruth, false);
    destroyPMat(QComp, false);

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}
