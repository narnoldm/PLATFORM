#include <dirent.h>

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

	int stdout_proc = 0;
	inputFile.getParamInt("stdout_proc", stdout_proc, 0);
	if (rank != stdout_proc) {
		cout.rdbuf(sink.rdbuf());
	}

	// ----- INPUT PARAMETERS -----

	string fieldInputString, samplingPointsFile, basisInputString, outDir;
    string centerFile, centerMethod, scaleFile, scaleMethod;
	bool center, normalize, outRegField, outAbsErrField;

	inputFile.getParamString("fieldInputString", fieldInputString);  // input token for FOM data series
	inputFile.getParamString("samplingPointsFile", samplingPointsFile);
	inputFile.getParamString("basisInputString", basisInputString);
	inputFile.getParamString("outDir", outDir);  // relative to field data directory

	// centering field data before regression
	inputFile.getParamBool("center", center);
	if (center) {
		// path to data centering profile
		// if not provided, use mean field
		inputFile.getParamString("centerFile", centerFile, "");
        inputFile.getParamString("centerMethod", centerMethod, "");
	}

	// normalizing field data before regression (after centering, if requested)
	inputFile.getParamBool("normalize", normalize);
	if (normalize) {
		// path to data normalization profile
		// if not provided, use normalization constants provided in fomInputString
		inputFile.getParamString("scaleFile", scaleFile, "");
        inputFile.getParamString("scaleMethod", scaleMethod, "");
	}

	inputFile.getParamBool("outRegField", outRegField, false);  // output unsteady regressed fields
	inputFile.getParamBool("outAbsErrField", outAbsErrField, false);  // output absolute unsteady error fields

	// ----- FINISH INPUT PARAMETERS -----

	// ----- PREPROCESSING -----

	PGrid *evenG;
	evenG = new PGrid(rank, size, 0);
	vector<string> token;

	tokenparse(fieldInputString, "|", token);
	tecIO* setField = new tecIO(token, "");

	token.clear();
	tokenparse(basisInputString, "|", token);
	meta* setBasis = new meta(token);

	size_t fieldDirPos = setField->prefix.find_last_of("/");
	string fieldDirBase = setField->prefix.substr(0, fieldDirPos);
	string fieldDir = fieldDirBase + "/regression/k" + to_string(setBasis->nSets) + "/" + outDir;
	if (rank == 0) {
		DIR* dir = opendir(fieldDirBase.c_str());
		// base directory does not exist
		if (!dir) {
			throw(-1);
		}
		size_t err = system(("mkdir -p " + fieldDir).c_str());
	}

	string errSuffix = "_" + to_string(setField->snap0) + "_" + to_string(setField->snapF) + "_" + to_string(setField->snapSkip);

	// ----- FINISH PREPROCESSING -----

	// ----- DATA LOADING -----

	// load field data
	pMat* QTruth = new pMat(setField->nPoints, setField->nSets, evenG, false);
	setField->batchRead(QTruth);

	// center data, if requested
	if (center) {
		if (centerFile == "") {
			setField->calcCentering(QTruth, centerMethod);
		} else {
			setField->calcCentering(QTruth, centerFile);
		}
		setField->centerData(QTruth);
	}

	// normalized data, if requested
	if (normalize) {
		if (scaleFile == "") {
			setField->calcScaling(QTruth, scaleMethod);
		} else {
			cout << "Norm file read not implemented yet" << endl;
			throw(-1);
		}
		setField->scaleData(QTruth);
	}

	// load regression basis
	pMat* basis = new pMat(setBasis->nPoints, setBasis->nSets, evenG, false);
	setBasis->batchRead(basis);

	// load sampling file
	cout << "Getting sample indices" << endl;
	fstream sampFile(samplingPointsFile, ios_base::in);
	int nSamps;
	sampFile >> nSamps;
	vector<int> sampleIdxs(nSamps);
	for (int i = 0; i < nSamps; ++i) {
		sampFile >> sampleIdxs[i];
		--sampleIdxs[i];
	}
	sampFile.close();

	// ----- FINISH DATA LOADING -----

	// ----- COMPUTE REGRESSION -----

	int nDOF = nSamps * setField->numVars;

	// extract rows of basis
	cout << "Extracting basis rows" << endl;
	pMat *basisSamp = new pMat(nDOF, setBasis->nSets, evenG, false);
	for (int i = 0; i < nSamps; i++) {
		cout << (double)i / nSamps * 100 << " percent points extracted \r";
		for (int j = 0; j < setField->numVars; j++)
			basisSamp->changeContext(basis, 1, setBasis->nSets, sampleIdxs[i] + j * setField->nCells, 0, i + j * nSamps, 0, false);
	}
	cout << endl;

	// compute pseudo-inverse of sampled basis, [S^T * U]^+
	cout << "Computing pseudo-inverse" << endl;
	pMat *pinvBasisSamp = new pMat(setBasis->nSets, nDOF, evenG, false);
	pinvBasisSamp->pinv(basisSamp);
	destroyPMat(basisSamp, false);

	// extract rows of field data, S^T * Q
	cout << "Extracting field rows" << endl;
	pMat *QTruthSamp = new pMat(nDOF, setField->nSets, evenG, false);
	for (int i = 0; i < nSamps; i++) {
		cout << (double)i / nSamps * 100 << " percent points extracted \r";
		for (int j = 0; j < setField->numVars; j++)
			QTruthSamp->changeContext(QTruth, 1, setField->nSets, sampleIdxs[i] + j * setField->nCells, 0, i + j * nSamps, 0, false);
	}
	cout << endl;

	// compute regression, U * [S^T * U]^+ * S^T * Q
	pMat* sampRegress = new pMat(setBasis->nSets, setField->nSets, evenG, false);
	sampRegress->matrix_Product('N', 'N', setBasis->nSets, setField->nSets, nDOF, pinvBasisSamp, 0, 0, QTruthSamp, 0, 0, 1.0, 0.0, 0, 0);
	destroyPMat(QTruthSamp, false);

	pMat* QComp = new pMat(setField->nPoints, setField->nSets, evenG, false);
	QComp->matrix_Product('N', 'N', setField->nPoints, setField->nSets, setBasis->nSets, basis, 0, 0, sampRegress, 0, 0, 1.0, 0.0, 0, 0);
	destroyPMat(sampRegress, false);
	destroyPMat(basis, false);

	// de-normalize and de-center, if normalization/centering was requested
	if (normalize) {
		setField->scaleData(QTruth, true);
		setField->scaleData(QComp, true);
	}
	if (center) {
		setField->centerData(QTruth, true);
		setField->scaleData(QComp, true);
	}

	// compute SVD of pseudo-inverse of sampled basis
	cout << "Computing SVD" << endl;
	int minDim = min(pinvBasisSamp->M, pinvBasisSamp->N);
	vector<double> S(minDim);
	pMat *U = new pMat(pinvBasisSamp->M, minDim, evenG, false);
	pMat *VT = new pMat(minDim, pinvBasisSamp->N, evenG, false);
	pinvBasisSamp->svd_run(pinvBasisSamp->M, pinvBasisSamp->N, 0, 0, U, VT, S, false);
	destroyPMat(U, false);
	destroyPMat(VT, false);
	destroyPMat(pinvBasisSamp, false);

	// ----- FINISH COMPUTE REGRESSION -----

	// ----- COMPUTE ERROR, WRITE OUTPUTS -----

	if (outRegField)
		setField->batchWrite(QComp, fieldDir, "reg_sol_");

	calc_abs_and_l2_error(QTruth, QComp, setField, fieldDir, errSuffix, outAbsErrField);

    // sampling error (maximum singular value)
	if (rank == 0) {
		cout << "Maximum singular value of [P^T U]^+: " << setprecision(16) << S[0] << endl;
		string outputFileName = fieldDir + "/samplingError.dat";
		ofstream outputFile(outputFileName);
		outputFile << setprecision(16) << S[0] << endl;
		outputFile.close();
	}

	// ----- FINISH COMPUTE ERROR, WRITE OUTPUTS -----

	destroyPMat(QTruth, false);
	destroyPMat(QComp, false);
	cout.rdbuf(strm_buffer);
	MPI_Finalize();
	return 0;
}
