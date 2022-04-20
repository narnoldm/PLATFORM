#include "metadata.hpp"
#include "param.hpp"

#include "sampling.hpp"

#include <unordered_set>

using namespace :: std;

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	ofstream sink("/dev/null");
	streambuf *strm_buffer = cout.rdbuf();

	paramMap inputFile("QR_pre.inp", rank); 	// input file

	// timing variables and output file
	double t0_start, t0_end, t1_start, t1_end, t2_start, t2_end;
	string timingOutput = "timings.dat";

	// clear timing file if it already exists
	ofstream out;
	out.open(timingOutput, ofstream::out | ofstream::trunc);
	out.close();

	t0_start = MPI_Wtime(); // full program timing

	int debug_proc = 0;
	inputFile.getParamInt("stdout_proc", debug_proc);
	if (rank != debug_proc)
    {
		cout.rdbuf(sink.rdbuf());
	}

	// Gappy POD regressor output format, based on GEMS ROM type
	// 0: Full residual ([P^T URes]^+, where URes is the full residual basis)
	// 1: Galerkin (USol^T URes * [P^T URes]^+, where URes is the RHS basis)
	// 2: Separate residual ([P^T URes]^+, [P^T USol]^+, and URes^T USol, where URes is the RHS basis)
	int regressorFormat;
	inputFile.getParamInt("regressorFormat", regressorFormat);

	// Residual/RHS basis path
	string inputRes;
	inputFile.getParamString("resInputString", inputRes); 	// PDP-demarcated input for residual basis

	// parse residual/RHS basis string
	vector<string> tokenRes;
	cout << "Residual input string is: " << inputRes << endl;
	tokenparse(inputRes, "|", tokenRes);

	int numModesRes;
	string firstFileBasis;
	// setup residual basis input, compare parameters against solution basis input
	meta *datasetRes;
	datasetRes = new meta(tokenRes);
	firstFileBasis = datasetRes->prefix + to_string(datasetRes->snap0) + datasetRes->suffix;
	numModesRes = datasetRes->nSets;

	int numModesMax;

    // field dimensions
    int nCells, nVars, nDOF;
	inputFile.getParamInt("nCells", nCells);
	inputFile.getParamInt("nVars", nVars);
	nDOF = nCells * nVars;

    PGrid* evenG;
	evenG = new PGrid(rank, size, 0);

	// repeat input reads for solution data if performing preprocessing for Galerkin ROM
	bool inputMatch = false;
	bool modesDiff = false;
	int numModesSol = 0;
	string firstFileSnaps;
	meta *datasetSol;
    tecIO *solScaling;
    tecIO *resScaling;
	if (regressorFormat > 0)
    {

		string inputSol;
		inputFile.getParamString("solInputString", inputSol); 	// PDP-demarcated input for solution data/basis

		vector<string> tokenSol;
		cout << "Solution input string is: " << inputSol << endl;
		tokenparse(inputSol, "|", tokenSol);

		// if input are identical, bases will just be copied from solution basis to RHS basis
		if (inputSol == inputRes)
        {
			cout << "Solution and RHS input strings were identical..." << endl;
			inputMatch = true;
		}

		// setup solution basis input
		datasetSol = new meta(tokenSol);
		firstFileSnaps = datasetSol->prefix + to_string(datasetSol->snap0) + datasetSol->suffix;
		numModesSol = datasetSol->nSets;

		int metaCheck = compareMeta(datasetSol, datasetRes);
		if (metaCheck == 2)
        {
			modesDiff = true;
			cout << "Solution and RHS bases have same datasets, but different mode counts..." << endl;
		}

		numModesMax = max(numModesSol, numModesRes); // determines which basis has more modes, only relevant if bases come from same dataset

        // read scaling data
        // TODO: this is SO JANK
        string scaleFile, scaleInput;
        bool scaleIsField;
        vector<string> token;

        // solution scaling
        scaleIsField = false;
        inputFile.getParamString("solScaleFile", scaleFile, "");
        inputFile.getParamString("solScaleInput", scaleInput, "");
        if (scaleFile != "")
        {
            scaleIsField = true;
        }
        else if (scaleInput != "")
        {
            scaleFile = "sphere";  // actual method is irrelevant, just has to be a valid method
        }
        else
        {
            cout << "Must provide either solScaleFile or solScaleInput is regressorFormat > 0" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        tokenparse(scaleInput, "|", token);
        solScaling = new tecIO(token);
        solScaling->calcScaling(NULL, scaleFile, scaleIsField, false);

        // residual/RHS scaling
        scaleIsField = false;
        inputFile.getParamString("resScaleFile", scaleFile, "");
        inputFile.getParamString("resScaleInput", scaleInput, "");
        if (scaleFile != "")
        {
            scaleIsField = true;
        }
        else if (scaleInput != "")
        {
            scaleFile = "sphere";  // actual method is irrelevant, just has to be a valid method
        }
        else
        {
            cout << "Must provide either resScaleFile or resScaleInput is regressorFormat > 0" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        tokenparse(scaleInput, "|", token);
        resScaling = new tecIO(token);
        resScaling->calcScaling(NULL, scaleFile, scaleIsField, false);

        // pre-compute P^-1 * G
        for (int i = 0; i < nDOF; ++i)
        {
            resScaling->scalingDivVec[i] = resScaling->scalingDivVec[i] / solScaling->scalingDivVec[i];
        }

	}
    else
    {
		numModesMax = numModesRes;
	}

	// percentage of total cells to be sampled
	double pSampling;
	inputFile.getParamDouble("pSampling", pSampling);

	pMat *USol, *USol_T;
	pMat *URes, *URes_T;

	// ##### SETTING UP RHS/RESIDUAL BASIS ##### //
	cout << "Loading RHS/residual basis..." << endl;

	// read modes from disk
	URes = new pMat(datasetRes->nPoints, datasetRes->nSets, evenG, false);

	t1_start = MPI_Wtime();
	if ( (regressorFormat > 0) && (inputMatch || (modesDiff && (numModesRes < numModesMax)))) {
		// if bases are identical, or if same basis set and solution basis has more modes, just load it now
		USol = new pMat(datasetSol->nPoints, datasetSol->nSets, evenG, false);
		datasetSol->batchRead(USol);
		URes->changeContext(USol, USol->M, numModesRes, 0, 0, 0, 0, false);
	} else {
		// if bases totally different, of solution basis set is same but has more modes than RHS basis, load solution basis now
		datasetRes->batchRead(URes);
	}
	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, timingOutput, "Residual/RHS basis load");

	// URes_T is just transpose of URes
	URes_T = new pMat(numModesRes, URes->M, evenG, false);
	URes_T->transpose(URes);

	// ##### FINISH RESIDUAL/RHS BASIS ##### //

	// ##### SETTING UP SOLUTION TRIAL BASIS ##### //
	if (regressorFormat > 0) {
		cout << "Loading solution basis..." << endl;

		// if inputMatch, basis has already been copied
		if (!inputMatch) {
			USol = new pMat(datasetSol->nPoints, datasetSol->nSets, evenG, false);

			// if same basis set, but solution basis has fewer modes than residual basis, extract first numModesSol modes from URes
			if (modesDiff && (numModesSol < numModesMax)) {
				USol->changeContext(URes, URes->M, numModesSol, 0, 0, 0, 0, false);

			// otherwise just load from disk
			} else {
				t1_start = MPI_Wtime();
				datasetSol->batchRead(USol);
				t1_end = MPI_Wtime();
				aggregateTiming(t1_end - t1_start, timingOutput, "Solution basis load");
			}
		}

		if (regressorFormat == 2) {
			USol_T = new pMat(numModesSol, USol->M, evenG, false);
			USol_T->transpose(USol);
		}

	}

	// ##### FINISH SOLUTION TRIAL BASIS ##### //

	// ##### START SAMPLING ##### //

	t1_start = MPI_Wtime();
	vector<int> gP; // will contain zero-indexed cell IDs of sampled cells
	unordered_set<int> samplingPoints;  // set version of gP, for automatically rejecting repeated entries
	int PointsNeeded = max(numModesMax, int(nCells * pSampling));
	int DOFNeeded = PointsNeeded * nVars;
    if (rank == 0) {
        cout << "Goal is " << PointsNeeded << " points" << endl;
    }

	// type of sampling
	// 0: QR only, no oversampling (QDEIM)
	// 1: QR + random oversampling (GappyPOD+R)
	// 2: QR + eigenvector-based oversampling (GappyPOD+E)
	// 3: DEIM greedy sampling (GappyPOD+D)
	int sampType;
	inputFile.getParamInt("sampType", sampType);

	// For greedy sampling methods, determine whether greedy objective function
	// is computed by cell or degree of freedom.
	int sampMethod;
	if (sampType > 1) {
		if (inputFile.getParamInt("sampMethod", sampMethod, 0)){
			switch(sampMethod) {
				case 0:
					cout << "Oversampling by DEGREE OF FREEDOM" << endl;
					break;
				case 1:
					cout << "Oversampling by CELL" << endl;
					break;
				default:
					cout << "Invalid value of sampMethod" << endl;
					throw(-1);
			}
		} else {
			cout << "Defaulting to oversampling by DEGREE OF FREEDOM" << endl;
		}
	}

    // ----- SEED SAMPLES ----- //

    string seedFile;
    if (inputFile.getParamString("seedFile", seedFile, "") && (rank == 0)) { 

		ifstream input(seedFile);
        string strIn;
        int numSeeds, idxIn;

        getline(input, strIn);
		numSeeds = atoi(strIn.c_str());

		// read and insert seed points
		cout << "Extracting " << numSeeds << " seed points..." << endl;
		for (int i = 0; i < numSeeds; i++) {

            getline(input, strIn);
            idxIn = atoi(strIn.c_str()) - 1;  // switch to 0 indexing

			cout << i << " " << idxIn << " " << endl;
			auto check = samplingPoints.emplace(idxIn);
			if (!check.second)
			{
				cout << "Repeated element" << endl;
			} else {
                gP.push_back(idxIn);
            }
		}
	}

	// ----- QR SAMPLING ----- //

	string qrSampBin_res, qrSampBin_sol;
	if ((sampType != 3) && (sampType != 4)) {
		t1_start = MPI_Wtime();

		qr_sampling(inputFile, "qrSampBinRes", "P_res", nCells, URes_T, gP, samplingPoints);
		if (regressorFormat == 2) {
			qr_sampling(inputFile, "qrSampBinSol", "P_sol", nCells, USol_T, gP, samplingPoints);
		}

		if (rank == 0) {
			cout << "Points after qr: " << samplingPoints.size() << " of " << PointsNeeded << endl;
		}
		MPI_Barrier(MPI_COMM_WORLD);
		t1_end = MPI_Wtime();
		aggregateTiming(t1_end - t1_start, timingOutput, "QR decomposition");
	}
	destroyPMat(URes_T, false);
	if (regressorFormat == 2)
		destroyPMat(USol_T, false);

	// ----- FINISH QR SAMPLING ----- //

	// ----- BOUNDARY SAMPLING ----- //

	// get boundary points
	t2_start = MPI_Wtime();

	// all sampling here is done by rank 0 process
	vector<int> itype;
	if (rank == 0) {

		// get desired boundary labels
		int numSampledBounds = 0;
		string labelInputString, percInputString;
		inputFile.getParamInt("numSampledBounds", numSampledBounds);
		vector<int> bcLabels(numSampledBounds);
		vector<double> bcPercs(numSampledBounds);
		for (int i = 0; i < numSampledBounds; ++i) {
			labelInputString = "boundLabel" + to_string(i+1);
			percInputString = "boundPerc" + to_string(i+1);
			inputFile.getParamInt(labelInputString, bcLabels[i]);
			cout << "Extracting boundary w/ label: " << bcLabels[i] << endl;
			if (inputFile.getParamDouble(percInputString, bcPercs[i], 1.0)) {
				cout << "Sampling " << to_string(100.0 * bcPercs[i]) << "% of boundary " << to_string(bcLabels[i]) << endl;
			} else {
				cout << "Assigning 100% sampling for boundary " << to_string(bcLabels[i]) << endl;
			}

		}

		// add boundary cells, if any
		if (numSampledBounds > 0) {
			cout << "Extracting boundary points..." << endl;
			// some I/O path strings
			string dfd_itype_file;
			inputFile.getParamString("dfd_itype_file", dfd_itype_file); // path to dfd_itype.bin generated by pgrid.x
			readMat(dfd_itype_file, itype);

			// count number of boundary cell types, get IDs of boundary cells
			vector<int> bcCount(numSampledBounds);
			vector<vector<int>> bcIDs(numSampledBounds);
			for (vector<int>::iterator it = itype.begin(); it != itype.end(); ++it) {
				for (int bc = 0; bc < numSampledBounds; ++bc) {
					if (bcLabels[bc] == *it) {
						bcCount[bc] += 1;
						bcIDs[bc].push_back(it - itype.begin());
						break;
					}
				}
			}

			// randomly sample boundary cells
			vector<int> bPoints;
			for (int bc = 0; bc < numSampledBounds; ++bc) {

				int numBoundsSamples = bcCount[bc] * bcPercs[bc];
				cout << "Sampling " << to_string(numBoundsSamples) << "/" << to_string(bcCount[bc]) << " cells from boundary " << to_string(bcLabels[bc]) << endl;
				bPoints.resize(bcCount[bc], 0);
				for (int i = 0; i < bPoints.size(); i++)
					bPoints[i] = i;
				srand(1); // seed random number generator
				random_shuffle(bPoints.begin(), bPoints.end());
				for (int it = 0; it < numBoundsSamples; ++it) {
					auto check = samplingPoints.emplace(bcIDs[bc][bPoints[it]]); //  add index of interator (i.e. cell_id, zero-indexed)
					if (!check.second) {
						cout << "Repeated element " << bcIDs[bc][bPoints[it]] << "\r";
					} else {
						gP.push_back(bcIDs[bc][bPoints[it]]); // need to keep track of this for GappyPOD+E
					}
				}
				cout << endl;
			}

			cout << endl << "Points after boundaries: " << samplingPoints.size() << " of " << PointsNeeded << endl;
		} else {
			cout << "No boundaries sampled..." << endl;
		}

        assert (PointsNeeded >= samplingPoints.size());
        cout << "Need " << PointsNeeded - samplingPoints.size() << " more points" << endl;

	}
	MPI_Barrier(MPI_COMM_WORLD);
	t2_end = MPI_Wtime();
	aggregateTiming(t2_end - t2_start, timingOutput, "Boundary sampling");

	// ----- FINISH BOUNDARY SAMPLING ----- //

	// ----- OVERSAMPLING ----- //

	t2_start = MPI_Wtime();

	// vector needed here to easily handle sampling based on multiple bases
	vector<pMat*> U_vec;
	U_vec.push_back(URes);
	if (regressorFormat == 2)
		U_vec.push_back(USol);

	switch (sampType) {

		// no oversampling
		case 0:
			cout << "No oversampling..." << endl;
			break;

		// random oversampling
		case 1:
			random_oversampling(nCells, PointsNeeded, samplingPoints, gP);
			break;

		// eigenvector-based oversampling
		case 2:
			eigenvector_oversampling(U_vec, sampMethod, nCells, nVars, PointsNeeded, samplingPoints, gP, timingOutput);
			break;

		// GNAT sampling, from Ben's preprint paper
		case 3:
			gnat_oversampling_peherstorfer(U_vec, sampMethod, nCells, nVars, PointsNeeded, samplingPoints, gP, timingOutput);
			break;

		// GNAT sampling, from Carlberg's 2017 paper
		case 4:
			gnat_oversampling_carlberg(U_vec, sampMethod, nCells, nVars, PointsNeeded, samplingPoints, gP, timingOutput);
			break;

		default:
			cout << "Invalid choice of sampling type: " << sampType << endl;
			throw(-1);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	t2_end = MPI_Wtime();
	aggregateTiming(t2_end - t2_start, timingOutput, "Oversampling");
	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, timingOutput, "Full sampling");

	// ----- FINISH OVERSAMPLING ----- //

	// write sampling points to disk
	if (rank == 0) {
		cout << "Writing sampling points to disk..." << endl;
		for_each(gP.begin(), gP.end(), [](int &tt) { tt += 1; }); 	// put in one-indexed format for writing to disk
		printASCIIVecP0("samplingPoints_unsorted.txt", gP, gP.size()); 		// writes unsorted sampling points to disk

        // sort gP and output ordered sampling points
        sort(gP.begin(), gP.end());
        printASCIIVecP0("samplingPoints.txt", gP, gP.size());
        writeMat("Pall.bin", gP.size(), 1, gP); 			// also writing as bin file
		for_each(gP.begin(), gP.end(), [](int &tt) { tt -= 1; }); 	// changing back to zero-indexed cell IDs
	}

	// communicate to all ranks
	if (rank != 0) {
		gP.resize(PointsNeeded, 0);
	}
	cout << "Broadcasting points..." << endl;
	MPI_Bcast(gP.data(), gP.size(), MPI_INT, 0, MPI_COMM_WORLD);

	// ##### FINISH SAMPLING ##### //

	// ##### OPERATOR PRECOMPUTING ##### //

	cout << "Calculating gappy POD regressors..." << endl;
	t1_start = MPI_Wtime();

	// [P^T * URes]^+
	pMat* pinvURes_samp = new pMat(numModesRes, DOFNeeded, evenG, false);
	calc_regressor(URes, pinvURes_samp, gP, nCells, nVars);

	pMat *regressorRes_out, *regressorSol_out, *USol_URes;

	if (regressorFormat > 0)
    {

        // [P^T * USol]^+
		if (regressorFormat == 2)
        {
			pMat* pinvUSol_samp = new pMat(numModesSol, DOFNeeded, evenG, false);
			calc_regressor(USol, pinvUSol_samp, gP, nCells, nVars);
            regressorSol_out = new pMat(DOFNeeded, numModesSol, evenG, false);
			regressorSol_out->transpose(pinvUSol_samp);
			destroyPMat(pinvUSol_samp, false);
		}

        // compute P^-1 *G * URes
        for (int i = 0; i < nDOF; ++i)
        {
            URes->scale_col_row(resScaling->scalingDivVec[i], i, true);
        }

        // USol^T * P^-1 * G * URes
		USol_URes = new pMat(numModesSol, numModesRes, evenG, false);
		USol_URes->matrix_Product('T', 'N', numModesSol, numModesRes, nDOF, USol, 0, 0, URes, 0, 0, 1.0, 0.0, 0, 0);

		// USol^T * P^-1 * G * URes * [P^T * URes]^+
		if (regressorFormat == 1) {

			pMat* projRegressor = new pMat(USol_URes->M, pinvURes_samp->N, evenG, false);
			projRegressor->matrix_Product('N', 'N', projRegressor->M, projRegressor->N, pinvURes_samp->M, USol_URes, 0, 0, pinvURes_samp, 0, 0, 1.0, 0.0, 0, 0);
			destroyPMat(USol_URes, false);
			destroyPMat(pinvURes_samp, false);

            regressorRes_out = new pMat(DOFNeeded, numModesSol, evenG, false);
			regressorRes_out->transpose(projRegressor);
			destroyPMat(projRegressor, false);
        }

	}

	// [P^T * URes]^+ for output
	if ((regressorFormat == 0) || (regressorFormat == 2)) {

        regressorRes_out = new pMat(DOFNeeded, numModesRes, evenG, false);
		regressorRes_out->transpose(pinvURes_samp);
		destroyPMat(pinvURes_samp, false);

	}

	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, "timings.dat", "Regressor calculation");

	destroyPMat(URes, false);
	if (regressorFormat > 0)
		destroyPMat(USol, false);

	// ##### FINISH OPERATOR PRECOMPUTING ##### //

	// ##### OUTPUT ##### //
	t1_start = MPI_Wtime();

	// write gappy POD regressors to disk

	meta* datasetResOut = new meta();
	if ((regressorFormat == 0) || (regressorFormat == 2)) {

		// write [P^T * URes]^+
		datasetResOut->snap0 = 1;
		datasetResOut->snapF = numModesRes;
		datasetResOut->snapSkip = 1;
		datasetResOut->nSets = numModesRes;
		datasetResOut->suffix = ".bin";
		datasetResOut->isInit = true;
		datasetResOut->nPoints = DOFNeeded;

		if (regressorFormat == 0) {

			// this represents the full residual regressor
			datasetResOut->batchWrite(regressorRes_out, "./", "pinv_res_");

		} else {

			// this represents the RHS regressor
			datasetResOut->batchWrite(regressorRes_out, "./", "pinv_rhs_");

			// write USol^T * URes
			USol_URes->write_bin("USolT_URes.bin");

			// write [P^T * USol]^+
			meta* datasetSolOut = new meta();
			datasetSolOut->snap0 = 1;
			datasetSolOut->snapF = numModesSol;
			datasetSolOut->snapSkip = 1;
			datasetSolOut->nSets = numModesSol;
			datasetSolOut->suffix = ".bin";
			datasetSolOut->isInit = true;
			datasetSolOut->nPoints = DOFNeeded;
			datasetSolOut->batchWrite(regressorSol_out, "./", "pinv_sol_");

		}

	} else {

		// write USol^T * URes * [P^T * URes]^+
		datasetResOut->snap0 = 1;
		datasetResOut->snapF = numModesSol;
		datasetResOut->snapSkip = 1;
		datasetResOut->nSets = numModesSol;
		datasetResOut->suffix = ".bin";
		datasetResOut->isInit = true;
		datasetResOut->nPoints = DOFNeeded;
		datasetResOut->batchWrite(regressorRes_out, "./", "pinv_rhs_");

	}

	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, "timings.dat", "Regressor mode write");

	// ##### FINISH OUTPUT ##### //

	t0_end = MPI_Wtime();
	aggregateTiming(t0_end - t0_start, timingOutput, "Complete run");

	// cleanup
	MPI_Barrier(MPI_COMM_WORLD);
	destroyPMat(regressorRes_out, false);
	if (regressorFormat == 2) {
		destroyPMat(regressorSol_out, false);
		destroyPMat(USol_URes, false);
	}

	cout.rdbuf(strm_buffer);
	MPI_Finalize();

	return 0;
}
