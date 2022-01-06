#include "metadata.hpp"
#include "param.hpp"

#include <set>

using namespace ::std;

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

	// DEIM interpolant output format, based on GEMS ROM type 
	// 0: SP-LSVT ([P^T U]^+)
	// 1: Galerkin (V^T U * [P^T U]^+) 
	int interpolantFormat;
	inputFile.getParamInt("interpolantFormat", interpolantFormat);

	// check whether user wants to compute bases or load them
	bool readRHSSnaps = false;
	inputFile.getParamBool("readRHSSnaps", readRHSSnaps); 	// whether to compute RHS basis here

	// RHS snapshots or basis paths
	string inputRHS;
	inputFile.getParamString("rhsInputString", inputRHS); 	// PDP-demarcated input for RHS data/basis

	// parse RHS basis string
	vector<string> tokenRHS;
	cout << "RHS input string is: " << inputRHS << endl;
	tokenparse(inputRHS, "|", tokenRHS);

	int numModesRHS;
	string firstFileBasis;
	// setup RHS/residual basis input, compare parameters against solution basis input
	meta *datasetRHS;
	if (readRHSSnaps) {
		datasetRHS = new tecIO(tokenRHS);
		firstFileBasis = datasetRHS->prefix + to_string(datasetRHS->snap0) + datasetRHS->suffix;
		dynamic_cast<tecIO *>(datasetRHS)->activateReorder(firstFileBasis.c_str());
		inputFile.getParamInt("numModesRHS", numModesRHS);

	// setup basis if reading basis
	} else {
		datasetRHS = new meta(tokenRHS);
		firstFileBasis = datasetRHS->prefix + to_string(datasetRHS->snap0) + datasetRHS->suffix;
		numModesRHS = datasetRHS->nSets;
	}

	int numModesMax;

	// repeat input reads for solution data if performing preprocessing for Galerkin ROM
	bool readSolSnaps = false;
	bool inputMatch = false;
	bool modesDiff = false;
	int numModesSol = 0;
	string firstFileSnaps;
	meta *datasetSol;
	if (interpolantFormat == 1) {
		
		inputFile.getParamBool("readSolSnaps", readSolSnaps); 	// whether to compute solution basis here

		string inputSol;
		inputFile.getParamString("solInputString", inputSol); 	// PDP-demarcated input for solution data/basis

		vector<string> tokenSol;
		cout << "Solution input string is: " << inputSol << endl;
		tokenparse(inputSol, "|", tokenSol);		

		// if input are identical, bases will just be copied from solution basis to RHS basis
		if (inputSol == inputRHS) {
			cout << "Solution and RHS input strings were identical..." << endl;
			inputMatch = true;
		}

		// setup solution basis input
		if (readSolSnaps) {
			datasetSol = new tecIO(tokenSol);
			firstFileSnaps = datasetSol->prefix + to_string(datasetSol->snap0) + datasetSol->suffix;
			dynamic_cast<tecIO *>(datasetSol)->activateReorder(firstFileSnaps.c_str());	
			inputFile.getParamInt("numModesSol", numModesSol); 	

			// if dataset is same and mode counts are not same, note this
			// I'm pretty sure the readRHSSnaps check is redundant, if inputMatch is already true
			if ( (readRHSSnaps) && (inputMatch) && (numModesSol != numModesRHS) ) {
				modesDiff = true;
				cout << "Solution and RHS bases have same datasets, but different mode counts..." << endl;
			}		 

		} else {
			datasetSol = new meta(tokenSol);
			firstFileSnaps = datasetSol->prefix + to_string(datasetSol->snap0) + datasetSol->suffix;
			numModesSol = datasetSol->nSets;

			if (!readRHSSnaps) {
				int metaCheck = compareMeta(datasetSol, datasetRHS);
				if (metaCheck == 2) {
					modesDiff = true;
					cout << "Solution and RHS bases have same datasets, but different mode counts..." << endl;
				} 
			}
		}

		numModesMax = max(numModesSol, numModesRHS); // determines which basis has more modes, only relevant if bases come from same dataset

	} else {

		numModesMax = numModesRHS;

	}

	// percentage of total cells to be sampled 
	double pSampling;
	inputFile.getParamDouble("pSampling", pSampling);

	// some I/O path strings
	string deimPrefix, dfd_itype_file;
	inputFile.getParamString("deimModesPrefix", deimPrefix); 	// prefix of DEIM interpolant mode files
	inputFile.getParamString("dfd_itype_file", dfd_itype_file); // path to dfd_itype.bin generated by pgrid.x

	PGrid *evenG;
	evenG = new PGrid(rank, size, 0);
	pMat *U, *A, *VT; 	// temporary basis matrices
	pMat *USol;
	pMat *URHS, *URHS_T;

	int nCells = 0, nVars = 0, nDOF = 0;

	// ##### SETTING UP RHS BASIS ##### //
	cout << "Loading RHS/residual basis..." << endl;

	// compute SVD of FOM snapshots, if computing POD basis here
	if (readRHSSnaps)
	{
		A = new pMat(datasetRHS->nPoints, datasetRHS->nSets, evenG, 0, 0, 0.0, false);
		t1_start = MPI_Wtime();
		datasetRHS->batchRead(A);
		t1_end = MPI_Wtime();
		aggregateTiming(t1_end - t1_start, timingOutput, "Residual snapshot load");

		nCells = dynamic_cast<tecIO *>(datasetRHS)->nCells;
		nVars  = dynamic_cast<tecIO *>(datasetRHS)->numVars;

		t1_start = MPI_Wtime();

		// TODO: this would need to include options to set centering and normalization profiles
		tecIO *datasetTec = dynamic_cast<tecIO *>(datasetRHS);
		datasetTec->calcAvg(A);
		datasetTec->subAvg(A);
		datasetTec->calcNorm(A);
		datasetTec->normalize(A);

		t1_end = MPI_Wtime();
		aggregateTiming(t1_end - t1_start, timingOutput, "Residual snapshot preprocessing");

		int MRHS = datasetRHS->nPoints; 
		int NRHS = datasetRHS->nSets;

		vector<double> SRHS;

		U = new pMat(MRHS, min(MRHS, NRHS), evenG, false);
		VT = new pMat(min(MRHS, NRHS), NRHS, evenG, false);
		SRHS.resize(min(MRHS, NRHS));

		A->svd_run(MRHS, NRHS, 0, 0, U, VT, SRHS);
		destroyPMat(A, false);
		destroyPMat(VT, false);

		// extract numModesRHS modes from U
		URHS = new pMat(U->M, numModesRHS, evenG, false);
		URHS->changeContext(U, U->M, numModesRHS, 0, 0, 0, 0, false);

		// if basis datasets match, also extract solution basis here
		if ( (interpolantFormat == 1) && inputMatch) {
			USol = new pMat(U->M, numModesSol, evenG, false);
			// if mode counts are different, extract modes from U
			if (modesDiff) {
				USol->changeContext(U, U->M, numModesSol, 0, 0, 0, 0, false);

			// if bases are identical
			} else {
				USol = URHS; // just use same address, don't need to allocate more memory
			}

		}

		destroyPMat(U, false); 
		
		// URHS_T is just transpose of URHS
		URHS_T = new pMat(numModesRHS, URHS->M, evenG, false);
		URHS_T->transpose(URHS, URHS_T->M, URHS_T->N, 0, 0);

	} else {

		// have to provide number of cells and variables, since reading from binary here (not SZPLT)
		inputFile.getParamInt("nCells", nCells); 
		inputFile.getParamInt("nVars", nVars);

		// read modes from disk
		URHS = new pMat(datasetRHS->nPoints, datasetRHS->nSets, evenG, false);

		t1_start = MPI_Wtime();
		if ( (interpolantFormat == 1) && (inputMatch || (modesDiff && (numModesRHS < numModesMax)))) {
			// if bases are identical, or if same basis set and solution basis has more modes, just load it now
			USol = new pMat(datasetSol->nPoints, datasetSol->nSets, evenG, false);
			datasetSol->batchRead(USol);
			URHS->changeContext(USol, USol->M, numModesRHS, 0, 0, 0, 0, false);
		} else {
			// if bases totally different, of solution basis set is same but has more modes than RHS basis, load solution basis now
			datasetRHS->batchRead(URHS);
		}
		t1_end = MPI_Wtime();
		aggregateTiming(t1_end - t1_start, timingOutput, "Residual mode load");

		// URHS_T is just transpose of URHS
		URHS_T = new pMat(numModesRHS, URHS->M, evenG, false); 						
		URHS_T->transpose(URHS);

	}

	nDOF = nCells * nVars;

	// ##### FINISH RHS BASIS ##### //

	// ##### SETTING UP SOLUTION TRIAL BASIS ##### // 
	if (interpolantFormat == 1) {
		cout << "Loading solution basis..." << endl;

		// compute SVD of FOM snapshots, if computing POD basis here
		if (readSolSnaps) {

			// if inputMatch, this has already been computed
			if (!inputMatch) {
				A = new pMat(datasetSol->nPoints, datasetSol->nSets, evenG, 0, 0, 0.0, false);
				t1_start = MPI_Wtime();
				datasetSol->batchRead(A);
				t1_end = MPI_Wtime();
				aggregateTiming(t1_end - t1_start, timingOutput, "Solution snapshot load");

				assert(nCells == dynamic_cast<tecIO *>(datasetSol)->nCells);
				assert(nVars  == dynamic_cast<tecIO *>(datasetSol)->numVars);

				t1_start = MPI_Wtime();

				// TODO: this would need to include options to set centering and normalization profiles
				tecIO *datasetTec = dynamic_cast<tecIO *>(datasetSol);
				datasetTec->calcAvg(A);
				datasetTec->subAvg(A);
				datasetTec->calcNorm(A);
				datasetTec->normalize(A);

				t1_end = MPI_Wtime();
				aggregateTiming(t1_end - t1_start, timingOutput, "Solution snapshot preprocessing");
				int MSol = datasetSol->nPoints; 
				int NSol = datasetSol->nSets;

				vector<double> SSol;

				U = new pMat(MSol, min(MSol, NSol), evenG, false);
				VT = new pMat(min(MSol, NSol), NSol, evenG, false);
				SSol.resize(min(MSol, NSol));

				A->svd_run(MSol, NSol, 0, 0, U, VT, SSol);
				destroyPMat(A, false);
				destroyPMat(VT, false);

				// extract numModesSol modes from U
				USol = new pMat(U->M, numModesSol, evenG, false);
				USol->changeContext(U, U->M, numModesSol, 0, 0, 0, 0, false);
				destroyPMat(U, false); // don't need this anymore

			}

		} else {

			// if inputMatch, basis has already been copied
			if (!inputMatch) {
				USol = new pMat(datasetSol->nPoints, datasetSol->nSets, evenG, false);

				// if same basis set, but solution basis has fewer modes than RHS basis, extract first numModesSol modes from URHS
				if (modesDiff && (numModesSol < numModesMax)) {
					USol->changeContext(URHS, URHS->M, numModesSol, 0, 0, 0, 0, false);

				// otherwise just load from disk
				} else {
					t1_start = MPI_Wtime();
					datasetSol->batchRead(USol);
					t1_end = MPI_Wtime();
					aggregateTiming(t1_end - t1_start, timingOutput, "Solution basis load");
				}
			}

		}
	}

	// ##### FINISH SOLUTION TRIAL BASIS ##### //

	// ##### START SAMPLING ##### //

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
	int greedyDOF;
	if (sampType > 1) {
		if(inputFile.getParamInt("sampMethod", sampMethod)){
			switch(sampMethod) {
				case 0:
					cout << "Oversampling by DEGREE OF FREEDOM" << endl;
					greedyDOF = nDOF;
					break;
				case 1:
					cout << "Oversampling by CELL" << endl;
					greedyDOF = nCells;
					break;
				default:
					cout << "Invalid value of sampMethod" << endl;
					throw(-1);
			}

		} else {
			sampMethod = 0;
			cout << "Defaulting to oversampling by DEGREE OF FREEDOM" << endl;
			greedyDOF = nDOF;
		}
	}

	// compute QR factorization of U^T
	// This writes the pivot indices to disk, since it's easier to do this than collect to rank 0 process
	int PointsNeeded = max(numModesRHS, int(nCells * pSampling)); // total number of cells that need to be sampled
	vector<int> P;
	string qrSampBin;
	if (sampType != 3) {
		if(inputFile.getParamString("qrSampBin", qrSampBin)){
			cout << "Retrieving QR sampling points from " << qrSampBin << endl;
		} else {
			cout << "qrSampBin not specified, computing it now." << endl;
			qrSampBin = "P.bin";
			t1_start = MPI_Wtime();
			cout << "Computing QR decomposition..." << endl;
			URHS_T->qr_run(URHS_T->M, URHS_T->N, 0, 0, P, "./", false);  // contents of URHS_T are DESTROYED during QR decomposition
			t1_end = MPI_Wtime();
			aggregateTiming(t1_end - t1_start, timingOutput, "QR decomposition");
		}
	}
	destroyPMat(URHS_T, false); 

	vector<int> gP; // will contain zero-indexed cell IDs of sampled cells
	vector<int> itype;
	set<int> samplingPoints;  // set version of gP, for automatically rejecting repeated entries

	// get QR and boundary points
	t1_start = MPI_Wtime();
	t2_start = MPI_Wtime();

	// all sampling here is done by rank 0 process
	if (rank == 0)
	{
		if (sampType != 3) {

			// read QR samples
			// gP MAY contain DOFs corresponding to the SAME CELL at this point
			readMat(qrSampBin, gP);
			gP.resize(numModesRHS);

			// sampled QR cells
			cout << "Extracting QR points..." << endl;
			for (int i = 0; i < gP.size(); i++)
			{
				gP[i]--; //switch to 0 indexing

				//switch to zero-indexed cell IDs
				cout << i << " " << gP[i] << " " << endl;
				gP[i] = gP[i] % nCells;
				auto check = samplingPoints.emplace(gP[i]);
				if (!check.second)
				{
					cout << "Repeated element" << endl;
				}
			}
			cout << "Goal is " << PointsNeeded << " points" << endl;
			cout << "Points after qr: " << samplingPoints.size() << " of " << PointsNeeded << endl;

			// gP now only contains unique sampled cell indices
			gP.resize(samplingPoints.size(), 0);
			copy(samplingPoints.begin(), samplingPoints.end(), gP.begin());

		}

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
			try {
				inputFile.getParamDouble(percInputString, bcPercs[i]);
				cout << "Sampling " << to_string(100.0 * bcPercs[i]) << "% of boundary " << to_string(bcLabels[i]) << endl;
			} catch (int e) {
				cout << "Assigning 100% sampling for boundary " << to_string(bcLabels[i]) << endl;
				bcPercs[i] = 1.0;
			}
			
		}

		// add boundary cells, if any
		if (numSampledBounds > 0) {
			cout << "Extracting boundary points..." << endl;
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
			}

			cout << endl << "Points after boundaries: " << samplingPoints.size() << " of " << PointsNeeded << endl;
			assert(PointsNeeded >= samplingPoints.size());
			cout << "Need " << PointsNeeded - samplingPoints.size() << " more points" << endl;
		} else {
			cout << "No boundaries sampled..." << endl;
		}

	}
	MPI_Barrier(MPI_COMM_WORLD);
	t2_end = MPI_Wtime(); // this is a serial timing, so has to be right after barrier
	aggregateTiming(t2_end - t2_start, timingOutput, "Boundary sampling");

	// oversampling
	bool allPoints = false;
	t2_start = MPI_Wtime();

	// random oversampling
	if (sampType == 0) {
		cout << "No oversampling..." << endl;
	} else if (sampType == 1) {
		// do everything on rank 0
		if (rank == 0) {

			cout << "Randomly oversampling..." << endl;

			srand(1);				// seed random number generator
			vector<int> rPoints;
			rPoints.resize(nCells, 0);
			for (int i = 0; i < rPoints.size(); i++)
				rPoints[i] = i;

			// selecting points from a randomly shuffled list of integers from 0 to ncells-1				
			random_shuffle(rPoints.begin(), rPoints.end());
			for (vector<int>::iterator it = rPoints.begin(); it != rPoints.end(); ++it)
			{
				auto check = samplingPoints.emplace(*it);
				if (!check.second)
				{
					cout << "repeated element " << *it << "\r";
				}
				if (samplingPoints.size() == PointsNeeded) {
					cout << endl << "All points found..." << endl;
					allPoints = true;
					break;
				}
			}

			// this should theoretically never happen
			if (!allPoints) {
				cout << endl << "Somehow, could not find enough sampled, this should never happen." << endl;
				return(-1);
			}

		}

	// eigenvector-based oversampling
	} else if (sampType == 2) {

		// declare some variables for this routine
		pMat *U_E, *VT_E;
		vector<double> S_E;
		int M, N, minDim;
		int numInitPoints, numCurrentPoints, numCurrentDOFs;
		int cellID;

		pMat *rVec = new pMat(1, nDOF, evenG, false);
		pMat *rVecCell;
		if (sampMethod == 1)
			rVecCell = new pMat(1, nCells, evenG, false);
		pMat *nonUniqueVec = new pMat(1, nDOF, evenG, false); // marker to determine if a degree of freedom has already been selected

		// sort and broadcast gP to all processes
		if (rank == 0)
			numInitPoints = gP.size();

		MPI_Bcast(&numInitPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank != 0)
			gP.resize(numInitPoints, 0);
		MPI_Bcast(gP.data(), gP.size(), MPI_INT, 0, MPI_COMM_WORLD);
		numCurrentPoints = numInitPoints;

		// mark DOFs that have already been sampled
		for (int j = 0; j < numInitPoints; ++j) {
			for (int k = 0; k < nVars; ++k) {
				nonUniqueVec->setElement(0, k * nCells + gP[j], 1.0);
			}
		}

		// Grab initial cell rows, BUT ordered by cell, then by var
		// This is done because newly sampled rows get appended to END of URHS_samp_E,
		// and we want to just copy the first N rows when copying URHS_samp_E
		pMat *URHS_samp_E = new pMat(PointsNeeded*nVars, datasetRHS->nSets, evenG, false);
		pMat *URHS_samp_E_copy = new pMat(PointsNeeded*nVars, datasetRHS->nSets, evenG, false);
		for (int j = 0; j < numCurrentPoints; ++j) {
			cout << (double)j / numCurrentPoints * 100 << " percent points extracted \r";
			for (int k = 0; k < nVars; ++k)
				URHS_samp_E_copy->changeContext(URHS, 1, numModesRHS, gP[j] + k * nCells, 0, j * nVars + k, 0, false);
		}
		cout << endl;

		// local timing variables
		double tE_start, tE_start2;
		double tE1, tE2, tE3, tE4, tE5, tE6, tE7, tE8, tE9, tE10;
		tE1 = 0.0; tE2 = 0.0; tE3 = 0.0; tE4 = 0.0; tE5 = 0.0; tE6 = 0.0;

		// variables for computing argsort
		int breakSignal;
		double maxLocal, maxGlobal;
		int argMaxLocal, argMaxGlobal;

		// loop over number of requires samples left
		int outFreq = (PointsNeeded - numInitPoints)/1000 + 1; 	// this is potentially a lot of output, scale to roughly 1000 lines of output
		for (int i = 0; i < (PointsNeeded - numInitPoints); ++i) {

			if ( (i % outFreq) == 0)
				cout << (double)i / (PointsNeeded - numInitPoints) * 100 << " percent GappyPOD+E points sampled \r";

			numCurrentDOFs = numCurrentPoints * nVars;

			// copy first numCurrentDOFs rows from URHS_samp_E_copy to URHS_samp_E
			// we do this because URHS_samp_E will be destroyed during the SVD
			tE_start = MPI_Wtime();
			URHS_samp_E->changeContext(URHS_samp_E_copy, numCurrentDOFs, numModesRHS, 0, 0, 0, 0, false);
			tE1 += (MPI_Wtime() - tE_start);

			// compute SVD of sampled URHS, hold on to singular values and RIGHT singular vectors
			M = numCurrentDOFs; 
			N = numModesRHS;
			minDim = min(M, N);

			U_E = new pMat(M, minDim, evenG, false);  // this will always change, must be reallocated/deleted
			VT_E = new pMat(minDim, N, evenG, false); // this *might* not need to be reallocated/deleted every iterations, could always be [numModesRHS x numModesRHS]
			S_E.resize(minDim, 0.0);

			tE_start = MPI_Wtime();
			URHS_samp_E->svd_run(M, N, 0, 0, U_E, VT_E, S_E, false); // contents of URHS_samp_E are destroyed during SVD
			tE2 += (MPI_Wtime() - tE_start);
			destroyPMat(U_E, false);

			// compute vector-matrix product of last right singular vector transposed and URHS transposed
			tE_start = MPI_Wtime();
			rVec->matrix_Product('N', 'T', 1, nDOF, minDim, VT_E, minDim-1, 0, URHS, 0, 0, 1.0, 0.0, 0, 0);
			tE3 += (MPI_Wtime() - tE_start);
			destroyPMat(VT_E, false);

			tE_start = MPI_Wtime();
			// zero out DOFs that have already been selected
			for (int j = 0; j < rVec->dataD.size(); ++j) {
				if (nonUniqueVec->dataD[j] == 1.0)
					rVec->dataD[j] = 0.0;
			}

			// square elements
			for (int j = 0; j < rVec->dataD.size(); ++j)
				rVec->dataD[j] = rVec->dataD[j] * rVec->dataD[j];

			// if sampling by cell, compute sum for each cell
			if (sampMethod == 1) {
				// zero out to start
				for (int j = 0; j < rVecCell->dataD.size(); ++j)
					rVecCell->dataD[j] = 0.0;
				for (int j = 0; j < nVars; ++j) {
					rVecCell->matrix_Sum('N', 1, nCells, rVec, 0, j * nCells, 0, 0, 1.0, 1.0);
				}
			}
			tE4 += (MPI_Wtime() - tE_start);

			// get unique cell ID
			tE_start = MPI_Wtime();

			breakSignal = 0;
			maxGlobal = 0;
			maxLocal = 0;
			argMaxLocal = -1;
			argMaxGlobal = -1;
			
			tE_start2 = MPI_Wtime();
			if (sampMethod == 1) {
				rVecCell->dMax(1, 0, maxLocal, argMaxLocal);
			} else {
				rVec->dMax(1, 0, maxLocal, argMaxLocal);
			}
			tE7 += (MPI_Wtime() - tE_start2);

			// Allreduce(MPI_MAX) the maximum value
			tE_start2 = MPI_Wtime();
			MPI_Allreduce(&maxLocal, &maxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
			tE8 += (MPI_Wtime() - tE_start2);

			// if rank doesn't own the max, set argMaxLocal = -1 so that MPI_MAX can determine the correct argMaxGlobal
			if (maxLocal != maxGlobal) {
				argMaxLocal = -1;
			}
			tE_start2 = MPI_Wtime();
			MPI_Allreduce(&argMaxLocal, &argMaxGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			tE9 += (MPI_Wtime() - tE_start2);

			// emplace cell ID in samplingPoints
			tE_start2 = MPI_Wtime();
			if (sampMethod == 1) {
				cellID = argMaxGlobal;
			} else {
				cellID = argMaxGlobal % nCells;
			}

			if (rank == 0) {
				auto check = samplingPoints.emplace(cellID);
				if (check.second) {
					breakSignal = 1;
				}
			}
			MPI_Allreduce(MPI_IN_PLACE, &breakSignal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			tE10 += (MPI_Wtime() - tE_start2);

			if (breakSignal == 0) {
				cout << "Non-unique cell found in greedy algorithm" << endl;
				cout << "Something went wrong..." << endl;
				cout << "argMaxGlobal: " << argMaxGlobal << endl;
				cout << "maxGlobal: " << maxGlobal << endl;
				MPI_Barrier(MPI_COMM_WORLD);
				throw(-1);
			}

			// mark all DOFs associated with selected cell
			for (int k = 0; k < nVars; ++k)
				nonUniqueVec->setElement(0, k * nCells + cellID, 1.0);

			tE5 += (MPI_Wtime() - tE_start);

			// insert into gP vector
			gP.push_back(cellID);
			numCurrentPoints++;

			// append newly sampled rows of URHS to URHS_samp_E_copy
			tE_start = MPI_Wtime();
			for (int k = 0; k < nVars; ++k)
				URHS_samp_E_copy->changeContext(URHS, 1, numModesRHS, k * nCells + gP.back(), 0, (numCurrentPoints-1) * nVars + k, 0, false);
			tE6 += (MPI_Wtime() - tE_start);

		}

		// aggregate timings
		aggregateTiming(tE1, timingOutput, "GPOD+E - URHS_samp_E copy");
		aggregateTiming(tE2, timingOutput, "GPOD+E - SVD");
		aggregateTiming(tE3, timingOutput, "GPOD+E - rVec solve");
		aggregateTiming(tE4, timingOutput, "GPOD+E - Zero and square");
		aggregateTiming(tE5, timingOutput, "GPOD+E - Find unique");
		aggregateTiming(tE6, timingOutput, "GPOD+E - Append rows");
		aggregateTiming(tE7, timingOutput, "GPOD+E - pdamax_");
		aggregateTiming(tE8, timingOutput, "GPOD+E - max allreduce");
		aggregateTiming(tE9, timingOutput, "GPOD+E - argmax reduce");
		aggregateTiming(tE10, timingOutput, "GPOD+E - break allreduce");

		// cleanup
		cout << endl;
		destroyPMat(URHS_samp_E, false);
		destroyPMat(URHS_samp_E_copy, false);
		destroyPMat(rVec, false);
		destroyPMat(nonUniqueVec, false);

	// DEIM greedy sampling
	} else if (sampType == 3) {

		// declare some variables for this routine
		int numCurrentDOFs, numInitPoints;
		int modeThresh, modeIdx;
		int cellID;
		
		pMat *rVec = new pMat(nDOF, 1, evenG, false);
		pMat *rVecCell;
		if (sampMethod == 1)
			rVecCell = new pMat(nCells, 1, evenG, false);
		pMat *nonUniqueVec = new pMat(nDOF, 1, evenG, false); // marker to determine if a degree of freedom has already been selected
		pMat *URHS_samp_E = new pMat(PointsNeeded * nVars, datasetRHS->nSets, evenG, false);
		pMat *URHS_samp_E_copy = new pMat(PointsNeeded * nVars, datasetRHS->nSets, evenG, false);

		// this is the same dimension as URHS_samp_E because PDGELS requires that URHS_samp_E and lsSol have same row block size
		// no way to force blocking sizes in pMat currently
		pMat *lsSol = new pMat(PointsNeeded * nVars, datasetRHS->nSets, evenG, false);

		// sort and broadcast gP to all processes
		if (rank == 0)
			numInitPoints = gP.size();

		MPI_Bcast(&numInitPoints, 1, MPI_INT, 0, MPI_COMM_WORLD);
		if (rank != 0)
			gP.resize(numInitPoints, 0);
		MPI_Bcast(gP.data(), gP.size(), MPI_INT, 0, MPI_COMM_WORLD);

		// mark DOFs that have already been sampled
		for (int j = 0; j < numInitPoints; ++j) {
			for (int k = 0; k < nVars; ++k) {
				nonUniqueVec->setElement(k * nCells + gP[j], 0, 1.0);
			}
		}

		// local timing variables
		double tD_start;
		double tD1, tD2, tD3, tD4, tD5, tD6;
		tD1 = 0.0; tD2 = 0.0; tD3 = 0.0; tD4 = 0.0; tD5 = 0.0; tD6 = 0.0;

		// variables for computing argsort
		int breakSignal;
		double maxLocal, maxGlobal;
		int argMaxLocal, argMaxGlobal;

		// extract first column of URHS to rVec
		rVec->changeContext(URHS, nDOF, 1, 0, 0, 0, 0, false);

		// loop over number of desired sampling points (index i)
		int outFreq = (PointsNeeded - numInitPoints)/1000 + 1;
		tD_start = MPI_Wtime();
		for (int i = 0; i < (PointsNeeded - numInitPoints); ++i) {

			if ( (i % outFreq) == 0)
				cout << (double)i / (PointsNeeded - numInitPoints) * 100 << " percent GappyPOD+D points sampled \r";

			tD_start = MPI_Wtime();
			// zero out DOFs that have already been selected
			for (int j = 0; j < rVec->dataD.size(); ++j) {
				if (nonUniqueVec->dataD[j] == 1.0)
					rVec->dataD[j] = 0.0;
			}

			// compute absolute value of rVec
			for (int j = 0; j < rVec->dataD.size(); ++j) {
				rVec->dataD[j] = abs(rVec->dataD[j]);
			}

			// if sampling by cell, compute sum for each cell
			if (sampMethod == 1) {
				// zero out to start
				for (int j = 0; j < rVecCell->dataD.size(); ++j)
					rVecCell->dataD[j] = 0.0;
				for (int j = 0; j < nVars; ++j) {
					rVecCell->matrix_Sum('N', nCells, 1, rVec, j * nCells, 0, 0, 0, 1.0, 1.0);
				}
			}
			tD1 += (MPI_Wtime() - tD_start);

			// find the argmax of rVec
			// loop until a unique point is added
			int breakSignal = 0;
			tD_start = MPI_Wtime();

			maxGlobal = 0;
			maxLocal = 0;
			argMaxLocal = -1;
			argMaxGlobal = -1;
			if (sampMethod == 1) {
				rVecCell->dMax(0, 0, maxLocal, argMaxLocal);
			} else {
				rVec->dMax(0, 0, maxLocal, argMaxLocal);
			}

			// Allreduce(MPI_MAX) the maximum value
			MPI_Allreduce(&maxLocal, &maxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			// if rank doesn't own the max, set argMaxLocal = -1 so that MPI_MAX can determine the correct argMaxGlobal
			if (maxLocal != maxGlobal) {
				argMaxLocal = -1;
			}
			MPI_Allreduce(&argMaxLocal, &argMaxGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			// emplace cell ID in samplingPoints, if unique entry found signal to all processes to break
			if (sampMethod == 1) {
				cellID = argMaxGlobal;
			} else {
				cellID = argMaxGlobal % nCells;
			}

			if (rank == 0) {
				auto check = samplingPoints.emplace(cellID);
				if (check.second) {
					breakSignal = 1;
				}
			}

			MPI_Allreduce(MPI_IN_PLACE, &breakSignal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
			tD2 += (MPI_Wtime() - tD_start);

			if (breakSignal == 0) {
				cout << "Non-unique cell found in greedy algorithm" << endl;
				cout << "Something went wrong..." << endl;
				cout << "argMaxGlobal: " << argMaxGlobal << endl;
				cout << "maxGlobal: " << maxGlobal << endl;
				MPI_Barrier(MPI_COMM_WORLD);
				throw(-1);
			}

			// mark all DOFs associated with selected cell
			for (int k = 0; k < nVars; ++k)
				nonUniqueVec->setElement(k * nCells + cellID, 0, 1.0);

			// insert into gP vector
			gP.push_back(cellID);
			numCurrentDOFs = (i+1) * nVars;

			// append newly sampled rows of URHS to URHS_samp_E_copy, and transfer to URHS_samp_E (will be destroyed in least-squares solve)
			tD_start = MPI_Wtime();
			for (int k = 0; k < nVars; ++k)
				URHS_samp_E_copy->changeContext(URHS, 1, numModesRHS, k * nCells + gP.back(), 0, i * nVars + k, 0, false);
			URHS_samp_E->changeContext(URHS_samp_E_copy, numCurrentDOFs, numModesRHS, 0, 0, 0, 0, false);
			tD3 += (MPI_Wtime() - tD_start);

			modeThresh = min(i+1, numModesRHS); // d in paper, forces ceiling of numModesRHS
			modeIdx = (i+1) % (numModesRHS-1); 	// k in paper, cycles through modes

			// set up inputs to least-squares solve
			lsSol->changeContext(URHS_samp_E, numCurrentDOFs, 1, 0, modeIdx, 0, 0, false);

			// solve least squares, should ALWAYS be an overdetermined system for systems with more than one variable
			// on exit, solution is in first modeThresh rows of lsSol
			// this is the c vector in paper
			tD_start = MPI_Wtime();
			lsSol->leastSquares('N', numCurrentDOFs, modeThresh, 1, URHS_samp_E, 0, 0, 0, 0);
			tD4 += (MPI_Wtime() - tD_start);

			// compute new rvec
			// first, U[:, 1:d]*c
			tD_start = MPI_Wtime();
			rVec->matrix_Product('N', 'N', nDOF, 1, modeThresh, URHS, 0, 0, lsSol, 0, 0, 1.0, 0.0, 0, 0);
			tD5 += (MPI_Wtime() - tD_start);

			// second, U[:, k] - U[:, 1:d]*c
			tD_start = MPI_Wtime();
			rVec->matrix_Sum('N', nDOF, 1, URHS, 0, modeIdx, 0, 0, 1.0, -1.0);
			tD6 += (MPI_Wtime() - tD_start);

		}

		// aggregate timings
		aggregateTiming(tD1, timingOutput, "GPOD+D - Zero and absolute val");
		aggregateTiming(tD2, timingOutput, "GPOD+D - Fine unique");
		aggregateTiming(tD3, timingOutput, "GPOD+D - Append rows");
		aggregateTiming(tD4, timingOutput, "GPOD+D - Least squares");
		aggregateTiming(tD5, timingOutput, "GPOD+D - rVec matmul");
		aggregateTiming(tD6, timingOutput, "GPOD+D - rVec matsum");

		// cleanup 
		cout << endl;
		destroyPMat(URHS_samp_E_copy, false);
		destroyPMat(URHS_samp_E, false);
		destroyPMat(rVec, false);
		destroyPMat(nonUniqueVec, false);
			
	} else {
		cout << "Invalid choice of sampling type: " << sampType << endl;
	} 
	MPI_Barrier(MPI_COMM_WORLD);
	t2_end = MPI_Wtime();
	aggregateTiming(t2_end - t2_start, timingOutput, "Oversampling");
	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, timingOutput, "Full sampling");

	if (rank == 0) {
		gP.resize(samplingPoints.size(), 0);
		vector<double> gPD; // gP, but doubles
		gPD.resize(nCells, 0.0);
		copy(samplingPoints.begin(), samplingPoints.end(), gP.begin());

		cout << "Writing sampling points to disk..." << endl;
		for_each(gP.begin(), gP.end(), [](int &tt) { tt += 1; }); 	// put in one-indexed format for writing to disk
		printASCIIVecP0("samplingPoints.txt", gP, gP.size()); 		// writes sampling points to disk
		writeMat("Pall.bin", gP.size(), 1, gP); 			// also writing as bin file
		for_each(gP.begin(), gP.end(), [](int &tt) { tt -= 1; }); 	// changing back to zero-indexed cell IDs

		// mark every sampled cell in gPD with 1.0
		for (int i = 0; i < samplingPoints.size(); i++)
			gPD[gP[i]] = 1.0;

		vector<string> Pname;
		string tempname = "sampling";
		Pname.push_back(tempname);

		// could probably generalize this without the conditionals
		if (readSolSnaps) {
			dynamic_cast<tecIO *>(datasetSol)->writeSingleFile("sampling.szplt", Pname, gPD.data(), firstFileSnaps);
		} else if (readRHSSnaps) {
			dynamic_cast<tecIO *>(datasetRHS)->writeSingleFile("sampling.szplt", Pname, gPD.data(), firstFileSnaps);
		}

	}

	if (rank != 0)
	{
		gP.resize(PointsNeeded, 0);
	}
	cout << "Broadcasting points..." << endl;
	MPI_Bcast(gP.data(), gP.size(), MPI_INT, 0, MPI_COMM_WORLD);

	cout << "Calculating DEIM interpolant..." << endl;
	
	// URHS_samp is P^T * URHS
	pMat *URHS_samp;
	URHS_samp = new pMat(gP.size() * nVars, numModesRHS, evenG, false);

	// extract rows of URHS
	// would it not be easier to just construct the full selection matrix P? Then multiply P^T * U?

	// permute rows with pdlapiv, then copy first gp.size()*nVars rows

	// multiply URHS by actual selection matrix

	// copying individual rows of URHS
	t1_start = MPI_Wtime();
	for (int i = 0; i < gP.size(); i++)
	{
		cout << (double)i / gP.size() * 100 << " percent points extracted \r";
		for (int j = 0; j < nVars; j++)
			URHS_samp->changeContext(URHS, 1, numModesRHS, gP[i] + j * nCells, 0, i + j * gP.size(), 0, false);
	}
	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, "timings.dat", "Residual basis row extraction");


	t1_start = MPI_Wtime();
	// compute [P^T * URHS]^+ component of DEIM interpolant
	cout << "Computing pseudo-inverse..." << endl;
	pMat *pinvURHS_samp = new pMat(URHS_samp->N, URHS_samp->M, evenG, false);
	pinvURHS_samp->pinv(URHS_samp);
	destroyPMat(URHS_samp, false);

	pMat *deimInterp_T;
	if (interpolantFormat == 1) {
		// compute USol^T * URHS
		cout << "Computing V^T * U ..." << endl;
		pMat *USol_URHS = new pMat(numModesSol, numModesRHS, evenG, false);
		USol_URHS->matrix_Product('T', 'N', numModesSol, numModesRHS, USol->M, USol, 0, 0, URHS, 0, 0, 1.0, 0.0, 0, 0); 

		// compute full DEIM interpolant for Galerkin, USol^T * URHS * [P^T * URHS]^+
		cout << "Computing complete DEIM interpolant..." << endl;
		pMat *deimInterp = new pMat(USol_URHS->M, pinvURHS_samp->N, evenG, false);
		deimInterp->matrix_Product('N', 'N', deimInterp->M, deimInterp->N, pinvURHS_samp->M, USol_URHS, 0, 0, pinvURHS_samp, 0, 0, 1.0, 0.0, 0, 0);
		

		destroyPMat(USol_URHS, false);
		destroyPMat(pinvURHS_samp, false);

		deimInterp_T = new pMat(deimInterp->N, deimInterp->M, deimInterp->pG, false);
		deimInterp_T->transpose(deimInterp);
		destroyPMat(deimInterp, false);

	} else {

		// full DEIM interpolant for SP-LSVT is just [P^T * URHS]^+
		deimInterp_T = new pMat(pinvURHS_samp->N, pinvURHS_samp->M, pinvURHS_samp->pG, false);
		deimInterp_T->transpose(pinvURHS_samp);
		destroyPMat(pinvURHS_samp, false);

	}
	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, "timings.dat", "Interpolant calculation");
	
	// write DEIM interpolant to disk
	// insert zeros where it is not sampled, for GEMS input
	cout << "Emplacing zeros into DEIM interpolant..." << endl;
	t1_start = MPI_Wtime();
	pMat *deimInterpOut = new pMat(datasetRHS->nPoints, numModesRHS, evenG, false);
	for (int i = 0; i < gP.size(); i++)
	{
		cout << (double)i / gP.size() * 100 << " percent points emplaced \r";
		for (int j = 0; j < nVars; j++)
			deimInterpOut->changeContext(deimInterp_T, 1, numModesRHS, i + j * gP.size(), 0, gP[i] + j * nCells, 0, false);
	}
	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, "timings.dat", "Interpolant zero emplacement");

	// what is the purpose of this? We can't use this stuff with GEMS, and it will get activated no matter what if computing the SVD
	t1_start = MPI_Wtime();
	if (readSolSnaps || readRHSSnaps) {
		tecIO *datasetOut = new tecIO();
		datasetOut->snap0 = 1;
		datasetOut->snapF = deimInterpOut->N;
		datasetOut->snapSkip = 1;
		datasetOut->nSets = deimInterpOut->N;
		//datasetOut->prefix = "deimInterp";
		//datasetOut->suffix = ".szplt";
		datasetOut->isInit = true;

		if (readSolSnaps) {
			datasetOut->meshFile = datasetSol->prefix + to_string(datasetSol->snap0) + datasetSol->suffix;
			datasetOut->varName = dynamic_cast<tecIO *>(datasetSol)->varName;
			datasetOut->varIndex = dynamic_cast<tecIO *>(datasetSol)->varIndex;
		} else {
			datasetOut->meshFile = datasetRHS->prefix + to_string(datasetRHS->snap0) + datasetRHS->suffix;
			datasetOut->varName = dynamic_cast<tecIO *>(datasetRHS)->varName;
			datasetOut->varIndex = dynamic_cast<tecIO *>(datasetRHS)->varIndex;
		}
		
		datasetOut->fixedMesh = true;
		datasetOut->getDimNodes();
		datasetOut->numVars = datasetOut->varName.size();
		datasetOut->nPoints = datasetOut->nCells * datasetOut->numVars;

		datasetOut->activateReorder(firstFileSnaps.c_str());
		datasetOut->activateGEMSbin(firstFileSnaps.c_str());
		datasetOut->batchWrite(deimInterpOut, "./", deimPrefix);
		
	} else {
		meta *datasetOut = new meta();
		datasetOut->snap0 = 1;
		datasetOut->snapF = deimInterpOut->N;
		datasetOut->snapSkip = 1;
		datasetOut->nSets = deimInterpOut->N;
		datasetOut->prefix = "deimInterp";
		datasetOut->suffix = ".bin";
		datasetOut->isInit = true;

		datasetOut->nPoints = nCells * nVars;
		assert(datasetOut->nPoints == nCells * nVars);
		datasetOut->batchWrite(deimInterpOut, "./", deimPrefix);
	}
	t1_end = MPI_Wtime();
	aggregateTiming(t1_end - t1_start, "timings.dat", "Interpolant mode write");

	t0_end = MPI_Wtime();
	aggregateTiming(t0_end - t0_start, timingOutput, "Complete run");

	cout.rdbuf(strm_buffer);
	MPI_Finalize();

	return 0;
}
