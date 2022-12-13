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
    if (rank != debug_proc)
    {
        cout.rdbuf(sink.rdbuf());
    }

    // ----- INPUT PARAMETERS -----

    string cellIDFile, fomInputString, romInputString, basisInputString;
    string centerFile, centerMethod, scaleFile, scaleMethod;
    bool calcMags, centerIsField, scaleIsField;
    bool center, scale, outProjField, outLatentCode, outAbsErrField;
    int errType;
    string magsInput;
    vector<string> magsToken;

    // 1: compute error between FOM and projected FOM solution
    // 2: compute error between FOM and ROM solution
    // 3: compute error between projected FOM and ROM solution
    inputFile.getParamInt("errType", errType);

    if (!((errType > 0) && (errType <= 4))) {
        cout << "Invalid errType: " << errType << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // calculating magnitudes
    inputFile.getParamBool("calcMags", calcMags, false);
    if (calcMags)
    {
        inputFile.getParamString("magsInput", magsInput);
        tokenparse(magsInput, "|", magsToken);
    }


    inputFile.getParamString("fomInputString", fomInputString);  // input token for FOM data series
    inputFile.getParamString("cellIDFile", cellIDFile, "");

    // input token for ROM data series
    if (errType > 1)
        inputFile.getParamString("romInputString", romInputString);

    // if FOM projection is required
    bool projFOM = false;
    if ((errType == 1) || (errType == 3))
        projFOM = true;

    if (projFOM)
    {
        // input token for spatial mode series
        inputFile.getParamString("basisInputString", basisInputString);

        // centering FOM data before projection
        inputFile.getParamBool("center", center);
        if (center)
        {
            // path to data centering profile
            // if not provided, use mean field
            inputFile.getParamString("centerFile", centerFile, "");
            inputFile.getParamString("centerMethod", centerMethod, "");
            inputFile.getParamBool("centerIsField", centerIsField, false);
        }

        // normalizing FOM data before projection (after centering, if requested)
        inputFile.getParamBool("scale", scale);
        if (scale)
        {
            // path to data normalization profile
            // if not provided, use normalization constants provided in fomInputString
            inputFile.getParamString("scaleFile", scaleFile, "");
            inputFile.getParamString("scaleMethod", scaleMethod, "");
            inputFile.getParamBool("scaleIsField", scaleIsField, false);
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
    tecIO* setFOM = new tecIO(token, cellIDFile);
    setFOM->activateReorder();

    // ROM dataset
    tecIO* setROM;
    if (errType > 1)
    {
        token.clear();
        tokenparse(romInputString, "|", token);
        setROM = new tecIO(token, cellIDFile);
        setROM->activateReorder();

        // error checking
        assert (setFOM->dim == setROM->dim);
        assert (setFOM->nCells == setROM->nCells);
        assert (setFOM->numVars == setROM->numVars);
        assert (setFOM->nSets == setROM->nSets);
    }

    // basis dataset
    meta* setBasis;
    if (projFOM)
    {
        token.clear();
        tokenparse(basisInputString, "|", token);
        setBasis = new meta(token);

        // error checking
        assert (setBasis->nPoints == setFOM->nPoints);
    }

    // output directories
    size_t fomDirPos = setFOM->prefix.find_last_of("/");
    fomDir = setFOM->prefix.substr(0, fomDirPos);
    if (projFOM)
    {
        projDir = fomDir + "/projection/k" + to_string(setBasis->nSets);
        if (outProjField || outLatentCode)
        {
            if (!rank)
                size_t ierr = system(("mkdir " + projDir).c_str());
        }
    }

    if (errType > 1)
    {
        size_t romDirPos = setROM->prefix.find_last_of("/");
        romDir = setROM->prefix.substr(0, romDirPos);
        errDir = romDir;

        if (errType == 2)
        {
            errSuffix = "_vs_raw";
        }
        else
        {
            errSuffix = "_vs_proj";
            errDir += "/projection/k" + to_string(setBasis->nSets);
        }

        if (calcMags)
            errSuffix += "_mag";

        errSuffix += "_" + to_string(setROM->snap0) + "_" + to_string(setROM->snapF) + "_" + to_string(setROM->snapSkip);
    }
    else
    {
        errDir = fomDir + "/projection/k" + to_string(setBasis->nSets);
        if (calcMags)
            errSuffix += "_mag";
        errSuffix += "_" + to_string(setFOM->snap0) + "_" + to_string(setFOM->snapF) + "_" + to_string(setFOM->snapSkip);
    }
    if (!rank)
        size_t ierr = system(("mkdir -p " + errDir).c_str());

    // ----- FINISH PREPROCESSING -----

    // ----- DATA LOADING -----

    // load FOM data
    pMat* QTruth = new pMat(setFOM->nPoints, setFOM->nSets, evenG, false);
    setFOM->batchRead(QTruth);

    pMat *QComp, *QTruth_proj, *latentCode, *basis;
    // project FOM data if requested
    if (projFOM)
    {

        // center data, if requested
        if (center)
        {
            if (centerFile == "")
            {
                setFOM->calcCentering(QTruth, centerMethod, centerIsField, false);
            }
            else
            {
                setFOM->calcCentering(QTruth, centerFile, true, false);
            }
            setFOM->centerData(QTruth, false);
        }

        // scale data, if requested
        if (scale)
        {
            if (scaleFile == "")
            {
                setFOM->calcScaling(QTruth, scaleMethod, scaleIsField, false);
            }
            else
            {
                setFOM->calcScaling(QTruth, scaleFile, true, false);
            }
            setFOM->scaleData(QTruth, false);
        }

        if (outProjField)
            setFOM->batchWrite_bin(QTruth, projDir, "fom_sol_raw_");

        // load trial basis
        pMat* basis = new pMat(setBasis->nPoints, setBasis->nSets, evenG, false);
        setBasis->batchRead(basis);

        // project data to low-dimension
        latentCode = new pMat(setBasis->nSets, setFOM->nSets, evenG, false);
        QTruth_proj = new pMat(setFOM->nPoints, setFOM->nSets, evenG, false);
        latentCode->matrix_Product('T', 'N', basis->N, QTruth->N, basis->M, basis, 0, 0, QTruth, 0, 0, 1.0, 0.0, 0, 0);
        QTruth_proj->matrix_Product('N', 'N', basis->M, latentCode->N, basis->N, basis, 0, 0, latentCode, 0, 0, 1.0, 0.0, 0, 0);

        // de-scale and de-center, if normalization/centering was requested
        if (scale)
        {
            setFOM->scaleData(QTruth_proj, true);
            if (errType == 1)
                setFOM->scaleData(QTruth, true);
        }
        if (center)
        {
            setFOM->centerData(QTruth_proj, true);
            if (errType == 1)
                setFOM->centerData(QTruth, true);
        }

        if (errType == 1)
        {
            QComp = QTruth_proj;
        }
        else
        {
            // don't need original FOM solution any more
            destroyPMat(QTruth, false);
            QTruth = QTruth_proj;
        }

        // output projected field or latent code history, if requested
        if (outProjField)
            setFOM->batchWrite_bin(QTruth_proj, projDir, "proj_sol_");
        if (outLatentCode)
            latentCode->write_bin(projDir + "/latentCode.bin");

        destroyPMat(latentCode, false);

    }

    if (errType > 1)
        QComp = new pMat(setROM->nPoints, setROM->nSets, evenG, false);

    // load ROM data
    if (errType > 1)
    {
        setROM->batchRead(QComp);
    }

    // ----- FINISH DATA LOADING -----

    // ----- COMPUTE ERROR, WRITE OUTPUTS -----

    calc_abs_and_l2_error(QTruth, QComp, setFOM, errDir, errSuffix, outAbsErrField, calcMags, magsToken);

    // ----- FINISH COMPUTE ERROR, WRITE OUTPUTS -----

    // cleanup
    MPI_Barrier(MPI_COMM_WORLD);
    destroyPMat(QTruth, false);
    destroyPMat(QComp, false);

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}
