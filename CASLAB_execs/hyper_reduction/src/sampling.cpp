#include "sampling.hpp"

using namespace :: std;

// computes sampling points from pivots of QR decomposition of U_T
void qr_sampling(paramMap inputFile, const string& qrSampFileStr, int nCells, pMat* U_T, vector<int>& gP, set<int>& samplingPoints) {

	string qrSampBin;
	if (inputFile.getParamString(qrSampFileStr, qrSampBin, "P.bin")){
			cout << "Retrieving QR sampling points from " << qrSampBin << endl;
	} else {
		cout << qrSampFileStr << " not specified, computing it now." << endl;
		vector<int> P;
		U_T->qr_run(U_T->M, U_T->N, 0, 0, P, "./", qrSampBin, false);  // contents of U_T are DESTROYED during QR decomposition
	}

	// read QR samples
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) {
		// gP MAY contain DOFs corresponding to the SAME CELL at this point
		readMat(qrSampBin, gP);
		gP.resize(U_T->M);

		// sampled QR cells
		cout << "Extracting QR points..." << endl;
		for (int i = 0; i < gP.size(); i++) {
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
	}

}

// randomly samples cells
void random_oversampling(int nCells, int PointsNeeded, set<int>& samplingPoints) {

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// do everything on rank 0
	if (rank == 0) {

		cout << "Randomly oversampling..." << endl;

		srand(1);				// seed random number generator
		vector<int> rPoints;
		rPoints.resize(nCells, 0);
		for (int i = 0; i < rPoints.size(); i++)
			rPoints[i] = i;

		// selecting points from a randomly shuffled list of integers from 0 to ncells-1
		bool allPoints = false;
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
			throw(-1);
		}

	}
}

// eigenvector-based sampling from Peherstorfer et al., 2018 (preprint)
void eigenvector_oversampling(const vector<pMat*> U_vec, int sampMethod, int nCells, int nVars,
							  int PointsNeeded, set<int>& samplingPoints, vector<int>& gP, string& timingOutput) {

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// declare some variables for this routine
	int nDOF = nCells * nVars;
	int numInitPoints, numCurrentPoints, numCurrentDOFs;
	int cellID;

	int numModesMax = 0;
	for (int i = 0; i < U_vec.size(); ++i) {
		numModesMax = max(numModesMax, U_vec[i]->N);
	}

	PGrid* evenG = U_vec[0]->pG; // reuse process grid

	// greedy sampling metrics
	pMat *rVec = new pMat(1, nDOF, evenG, false);  // to be reused for each basis
	pMat *rVecSum = new pMat(1, nDOF, evenG, false);  // sum for all bases
	pMat *rVecCell;  // metric for mesh cell
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
	// This is done because newly sampled rows get appended to END of U_samp,
	// and we want to just copy the first N rows when copying U_samp
	vector<pMat*> U_samp_vec;
	vector<pMat*> U_samp_copy_vec;
	for (int i = 0; i < U_vec.size(); ++i) {
		U_samp_vec.push_back(new pMat(PointsNeeded*nVars, U_vec[i]->N, evenG, false));
		U_samp_copy_vec.push_back(new pMat(PointsNeeded*nVars, U_vec[i]->N, evenG, false));
		for (int j = 0; j < numCurrentPoints; ++j) {
			cout << (double)j / numCurrentPoints * 100 << " percent points extracted \r";
			for (int k = 0; k < nVars; ++k)
				U_samp_copy_vec[i]->changeContext(U_vec[i], 1, U_vec[i]->N, gP[j] + k * nCells, 0, j * nVars + k, 0, false);
		}
		cout << endl;
	}

	// local timing variables
	double tE_start;
	double tE1, tE2, tE3;
	tE1 = 0.0; tE2 = 0.0; tE3 = 0.0;

	// variables for computing argsort
	int argMaxGlobal;

	// loop over number of requires samples left
	int outFreq = (PointsNeeded - numInitPoints)/1000 + 1; 	// this is potentially a lot of output, scale to roughly 1000 lines of output
	for (int i = 0; i < (PointsNeeded - numInitPoints); ++i) {

		if ( (i % outFreq) == 0)
			cout << (double)i / (PointsNeeded - numInitPoints) * 100 << " percent GappyPOD+E points sampled \r";

		numCurrentDOFs = numCurrentPoints * nVars;

		// zero out metric vector
		for (int j = 0; j < rVecSum->dataD.size(); ++j)
			rVecSum->dataD[j] = 0.0;

		// compute greedy sampling metric
		tE_start = MPI_Wtime();
		for (int j = 0; j < U_vec.size(); ++j) {
			eigenvector_oversampling_metric(U_vec[j], U_samp_vec[j], U_samp_copy_vec[j], rVec, rVecSum, nonUniqueVec, numCurrentDOFs);
		}

		// if sampling by cell, compute sum for each cell
		if (sampMethod == 1) {
			// zero out to start
			for (int j = 0; j < rVecCell->dataD.size(); ++j)
				rVecCell->dataD[j] = 0.0;
			for (int j = 0; j < nVars; ++j) {
				rVecCell->matrix_Sum('N', 1, nCells, rVecSum, 0, j * nCells, 0, 0, 1.0, 1.0);
			}
		}
		tE1 += (MPI_Wtime() - tE_start);

		// get unique cell ID
		tE_start = MPI_Wtime();
		if (sampMethod == 1) {
			argMaxGlobal = rVecCell->argmax_vec();
		} else {
			argMaxGlobal = rVecSum->argmax_vec();
		}

		// emplace cell ID in samplingPoints
		if (sampMethod == 1) {
			cellID = argMaxGlobal;
		} else {
			cellID = argMaxGlobal % nCells;
		}

		// check for uniqueness
		if (rank == 0) {
			auto check = samplingPoints.emplace(cellID);
			if (!check.second) {
				cout << "Non-unique cell found in greedy algorithm" << endl;
				cout << "Something went wrong..." << endl;
				cout << "argMaxGlobal: " << argMaxGlobal << endl;
				throw(-1);
			}
		}

		// mark all DOFs associated with selected cell
		for (int k = 0; k < nVars; ++k)
			nonUniqueVec->setElement(0, k * nCells + cellID, 1.0);

		tE2 += (MPI_Wtime() - tE_start);

		// insert into gP vector
		gP.push_back(cellID);
		numCurrentPoints++;

		// append newly sampled rows of U to U_samp_copy
		tE_start = MPI_Wtime();
		for (int j = 0; j < U_vec.size(); ++j) {
			for (int k = 0; k < nVars; ++k)
				U_samp_copy_vec[j]->changeContext(U_vec[j], 1, U_vec[j]->N, k * nCells + gP.back(), 0, (numCurrentPoints-1) * nVars + k, 0, false);
		}
		tE3 += (MPI_Wtime() - tE_start);

	}

	// aggregate timings
	aggregateTiming(tE1, timingOutput, "GPOD+E - Metric calculation");
	aggregateTiming(tE2, timingOutput, "GPOD+E - Find unique");
	aggregateTiming(tE3, timingOutput, "GPOD+E - Append rows");

	// cleanup
	cout << endl;
	for (int i = 0; i < U_vec.size(); ++i) {
		destroyPMat(U_samp_vec[i], false);
		destroyPMat(U_samp_copy_vec[i], false);
	}
	destroyPMat(rVec, false);
	destroyPMat(rVecSum, false);
	if (sampMethod == 1)
		destroyPMat(rVecCell, false);
	destroyPMat(nonUniqueVec, false);
}

// computes greedy sampling metric for a single basis via eigenvector-based oversampling
void eigenvector_oversampling_metric(pMat* U, pMat* U_samp, pMat* U_samp_copy, pMat* rVec, pMat* rVecSum, pMat* nonUniqueVec, int numCurrentDOFs) {

	// copy first numCurrentDOFs rows from URHS_samp_E_copy to URHS_samp_E
	// we do this because URHS_samp_E will be destroyed during the SVD
	U_samp->changeContext(U_samp_copy, numCurrentDOFs, U->N, 0, 0, 0, 0, false);

	// compute SVD of sampled URHS, hold on to singular values and RIGHT singular vectors
	int M = numCurrentDOFs;
	int N = U->N;
	int minDim = min(M, N);

	pMat* U_E = new pMat(M, minDim, U->pG, false);  // this will always change, must be reallocated/deleted
	pMat* VT_E = new pMat(minDim, N, U->pG, false); // this *might* not need to be reallocated/deleted every iterations
	vector<double> S_E(minDim, 0.0);

	U_samp->svd_run(M, N, 0, 0, U_E, VT_E, S_E, false); // contents of U_samp are destroyed during SVD
	destroyPMat(U_E, false);

	// compute vector-matrix product of last right singular vector transposed and U transposed
	rVec->matrix_Product('N', 'T', 1, U->M, minDim, VT_E, minDim-1, 0, U, 0, 0, 1.0, 0.0, 0, 0);
	destroyPMat(VT_E, false);

	for (int j = 0; j < rVec->dataD.size(); ++j) {
		if (nonUniqueVec->dataD[j] != 1.0) {
			// add contribution from sampled degrees of freedom
			rVecSum->dataD[j] += rVec->dataD[j] * rVec->dataD[j];
		}
	}

}

// GNAT sampling based on algorithm from Peherstorfer et al., 2018 (preprint)
void gnat_oversampling_peherstorfer(pMat* URes, pMat* USol, int sampMethod, int nCells, int nVars, int nDOF, int numModesRHS,
							  int PointsNeeded, set<int>& samplingPoints, vector<int>& gP, string& timingOutput) {

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	PGrid* evenG = URes->pG; // reuse process grid

	// declare some variables for this routine
	int numCurrentDOFs, numInitPoints;
	int modeThresh, modeIdx;
	int cellID;

	pMat *rVec = new pMat(nDOF, 1, evenG, false);
	pMat *rVecCell;
	if (sampMethod == 1)
		rVecCell = new pMat(nCells, 1, evenG, false);
	pMat *nonUniqueVec = new pMat(nDOF, 1, evenG, false); // marker to determine if a degree of freedom has already been selected
	pMat *URHS_samp_E = new pMat(PointsNeeded * nVars, numModesRHS, evenG, false);
	pMat *URHS_samp_E_copy = new pMat(PointsNeeded * nVars, numModesRHS, evenG, false);

	// this is the same dimension as URHS_samp_E because PDGELS requires that URHS_samp_E and lsSol have same row block size
	// no way to force blocking sizes in pMat currently
	pMat *lsSol = new pMat(PointsNeeded * nVars, numModesRHS, evenG, false);

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
	int argMaxGlobal;

	// extract first column of URHS to rVec
	rVec->changeContext(URes, nDOF, 1, 0, 0, 0, 0, false);

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
		tD_start = MPI_Wtime();
		if (sampMethod == 1) {
			argMaxGlobal = rVecCell->argmax_vec();
		} else {
			argMaxGlobal = rVec->argmax_vec();
		}

		// emplace cell ID in samplingPoints, if unique entry found signal to all processes to break
		if (sampMethod == 1) {
			cellID = argMaxGlobal;
		} else {
			cellID = argMaxGlobal % nCells;
		}

		if (rank == 0) {
			auto check = samplingPoints.emplace(cellID);
			if (!check.second) {
				cout << "Non-unique cell found in greedy algorithm" << endl;
				cout << "Something went wrong..." << endl;
				cout << "argMaxGlobal: " << argMaxGlobal << endl;
				throw(-1);
			}
		}
		tD2 += (MPI_Wtime() - tD_start);

		// mark all DOFs associated with selected cell
		for (int k = 0; k < nVars; ++k)
			nonUniqueVec->setElement(k * nCells + cellID, 0, 1.0);

		// insert into gP vector
		gP.push_back(cellID);
		numCurrentDOFs = (i+1) * nVars;

		// append newly sampled rows of URHS to URHS_samp_E_copy, and transfer to URHS_samp_E (will be destroyed in least-squares solve)
		tD_start = MPI_Wtime();
		for (int k = 0; k < nVars; ++k)
			URHS_samp_E_copy->changeContext(URes, 1, numModesRHS, k * nCells + gP.back(), 0, i * nVars + k, 0, false);
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
		rVec->matrix_Product('N', 'N', nDOF, 1, modeThresh, URes, 0, 0, lsSol, 0, 0, 1.0, 0.0, 0, 0);
		tD5 += (MPI_Wtime() - tD_start);

		// second, U[:, k] - U[:, 1:d]*c
		tD_start = MPI_Wtime();
		rVec->matrix_Sum('N', nDOF, 1, URes, 0, modeIdx, 0, 0, 1.0, -1.0);
		tD6 += (MPI_Wtime() - tD_start);

	}

	// aggregate timings
	aggregateTiming(tD1, timingOutput, "GNAT (Peherstorfer) - Zero and absolute val");
	aggregateTiming(tD2, timingOutput, "GNAT (Peherstorfer) - Find unique");
	aggregateTiming(tD3, timingOutput, "GNAT (Peherstorfer) - Append rows");
	aggregateTiming(tD4, timingOutput, "GNAT (Peherstorfer) - Least squares");
	aggregateTiming(tD5, timingOutput, "GNAT (Peherstorfer) - rVec matmul");
	aggregateTiming(tD6, timingOutput, "GNAT (Peherstorfer) - rVec matsum");

	// cleanup
	cout << endl;
	destroyPMat(URHS_samp_E_copy, false);
	destroyPMat(URHS_samp_E, false);
	destroyPMat(rVec, false);
	destroyPMat(nonUniqueVec, false);
}

// GNAT sampling based on algorithm from Carlberg et al., 2017
void gnat_oversampling_carlberg(pMat* URes, pMat* USol, int sampMethod, int nCells, int nVars, int nDOF, int numModesRHS,
							  int PointsNeeded, set<int>& samplingPoints, vector<int>& gP, string& timingOutput) {

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	PGrid* evenG = URes->pG; // reuse process grid

	// declare some variables for this routine
	int numCurrentDOFs, numInitPoints;
	int cellID;

	pMat *rVec = new pMat(nDOF, 1, evenG, false);
	pMat *rVecCell;
	if (sampMethod == 1)
		rVecCell = new pMat(nCells, 1, evenG, false);
	pMat *nonUniqueVec = new pMat(nDOF, 1, evenG, false); // marker to determine if a degree of freedom has already been selected
	pMat *URHS_samp_E = new pMat(PointsNeeded * nVars, numModesRHS, evenG, false);
	pMat *URHS_samp_E_copy = new pMat(PointsNeeded * nVars, numModesRHS, evenG, false);

	// this is the same dimension as URHS_samp_E because PDGELS requires that URHS_samp_E and lsSol have same row block size
	// no way to force blocking sizes in pMat currently
	pMat *lsSol = new pMat(PointsNeeded * nVars, numModesRHS, evenG, false);

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
	int argMaxGlobal;

	// extract first column of URHS to rVec
	rVec->changeContext(URes, nDOF, 1, 0, 0, 0, 0, false);

	// carbon copy from paper
	int nc = numModesRHS;  // number working columns
	int na = PointsNeeded - numInitPoints;  // remaining points needed
	int nb = 0;  // counter
	int nit = min(nc, na);  // number of greedy iterations
	int nai_min = floor((double)na / nc);
	int nai;

	// greedy loop
	tD_start = MPI_Wtime();
	for (int i = 0; i < nit; ++i) {

		cout << (double)i / nit * 100 << " percent GNAT points sampled \r";

		tD_start = MPI_Wtime();

		for (int j = 0; j < rVec->dataD.size(); ++j) {
			// zero out DOFs that have already been selected
			if (nonUniqueVec->dataD[j] == 1.0) {
				rVec->dataD[j] = 0.0;
			// compute elementwise square of rVec
			} else {
				rVec->dataD[j] = rVec->dataD[j] * rVec->dataD[j];
			}
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

		// sample cell loop
		// need an extra sampling point for this loop
		nai = nai_min;
		if ((i + 1) <= (na % nc))
			nai++;

		for (int j = 0; j < nai; ++j) {

			tD_start = MPI_Wtime();
			if (sampMethod == 1) {
				argMaxGlobal = rVecCell->argmax_vec();
			} else {
				argMaxGlobal = rVec->argmax_vec();
			}

			// emplace cell ID in samplingPoints
			if (sampMethod == 1) {
				cellID = argMaxGlobal;
			} else {
				cellID = argMaxGlobal % nCells;
			}

			if (rank == 0) {
				auto check = samplingPoints.emplace(cellID);
				if (!check.second) {
					cout << "Non-unique cell found in greedy algorithm" << endl;
					cout << "Something went wrong..." << endl;
					cout << "argMaxGlobal: " << argMaxGlobal << endl;
					throw(-1);
				}
			}
			tD2 += (MPI_Wtime() - tD_start);

			// mark all DOFs associated with selected cell
			for (int k = 0; k < nVars; ++k) {
				nonUniqueVec->setElement(k * nCells + cellID, 0, 1.0);
				if (sampMethod == 0) {
					rVec->setElement(k * nCells + cellID, 0, 0.0);
				}
			}
			if (sampMethod == 1)
				rVecCell->setElement(cellID, 0, 0.0);

			// insert into gP vector
			gP.push_back(cellID);
			numCurrentDOFs = (i + 1) * nVars;

			// append newly sampled rows of URHS to URHS_samp_E_copy, and transfer to URHS_samp_E (will be destroyed in least-squares solve)
			tD_start = MPI_Wtime();
			for (int k = 0; k < nVars; ++k)
				URHS_samp_E_copy->changeContext(URes, 1, numModesRHS, k * nCells + gP.back(), 0, i * nVars + k, 0, false);
			URHS_samp_E->changeContext(URHS_samp_E_copy, numCurrentDOFs, numModesRHS, 0, 0, 0, 0, false);
			tD3 += (MPI_Wtime() - tD_start);

		}

		// have gotten all the samples we need
		if (i == (nit - 1))
			break;

		// set up inputs to least-squares solve
		// this is the b term in || Ax - b ||
		lsSol->changeContext(URHS_samp_E, numCurrentDOFs, 1, 0, i + 1, 0, 0, false);

		// solve least squares, should ALWAYS be an overdetermined system for systems with more than one variable
		// on exit, solution is in first modeThresh rows of lsSol
		// this is the alpha vector in paper
		tD_start = MPI_Wtime();
		lsSol->leastSquares('N', numCurrentDOFs, i + 1, 1, URHS_samp_E, 0, 0, 0, 0);
		tD4 += (MPI_Wtime() - tD_start);

		// compute new rvec
		// first, U[:, 1:nb] * alpha
		tD_start = MPI_Wtime();
		rVec->matrix_Product('N', 'N', nDOF, 1, i + 1, URes, 0, 0, lsSol, 0, 0, 1.0, 0.0, 0, 0);
		tD5 += (MPI_Wtime() - tD_start);

		// second, U[:, k] - U[:, 1:nb] * alpha
		tD_start = MPI_Wtime();
		rVec->matrix_Sum('N', nDOF, 1, URes, 0, i + 1, 0, 0, 1.0, -1.0);
		tD6 += (MPI_Wtime() - tD_start);

	}

	// aggregate timings
	aggregateTiming(tD1, timingOutput, "GNAT (Carlberg) - Zero and absolute val");
	aggregateTiming(tD2, timingOutput, "GNAT (Carlberg) - Find unique");
	aggregateTiming(tD3, timingOutput, "GNAT (Carlberg) - Append rows");
	aggregateTiming(tD4, timingOutput, "GNAT (Carlberg) - Least squares");
	aggregateTiming(tD5, timingOutput, "GNAT (Carlberg) - rVec matmul");
	aggregateTiming(tD6, timingOutput, "GNAT (Carlberg) - rVec matsum");

	// cleanup
	cout << endl;
	destroyPMat(URHS_samp_E_copy, false);
	destroyPMat(URHS_samp_E, false);
	destroyPMat(rVec, false);
	destroyPMat(nonUniqueVec, false);
}

// computes gappy POD regressor [P^T * U]^+ for saving to disk
void calc_regressor(pMat* U, pMat* regressor, vector<int>& gP, int nCells, int nVars) {

	pMat *U_samp;
	U_samp = new pMat(gP.size() * nVars, U->N, U->pG, false);

	assert (regressor->M == U_samp->N);
	assert (regressor->N == U_samp->M);

	// copying rows of U
	for (int i = 0; i < gP.size(); i++) {
		cout << (double)i / gP.size() * 100 << " percent points extracted \r";
		for (int j = 0; j < nVars; j++)
			U_samp->changeContext(U, 1, U->N, gP[i] + j * nCells, 0, i + j * gP.size(), 0, false);
	}
	cout << endl;

	// compute [P^T * U]^+
	cout << "Computing pseudo-inverse..." << endl;
	regressor->pinv(U_samp);
	destroyPMat(U_samp, false);

}

/**
 *  Put regressor into GEMS format, with zeros at unsampled degrees of freedom
 *
 *  If transA = 'T', AIn is assumed to be in its original [P^T * U]^+ format, and must be transposed first.
 *  AOut must already be allocated, initialized with zeros, and have the same dimension as the basis (U) it is derived from.
 *  Does not destroy AIn, that is the user's responsibility.
 */
void emplace_zeros(const char transA, pMat* AIn, pMat* AOut, vector<int>& gP, int nCells, int nVars) {

	pMat* regressor;
	int numModes;

	// check dimensions of AIn and AOut, set regressor
	if (transA == 'T') {
		assert (AIn->N == (gP.size() * nVars));
		numModes = AIn->M;
		regressor = new pMat(AIn->N, AIn->M, AIn->pG, false);
		regressor->transpose(AIn);
	} else if (transA == 'F') {
		assert(AIn->M == (gP.size() * nVars));
		numModes = AIn->N;
		regressor = AIn;
	} else {
		cout << "Invalid value of transA: " << transA << endl;
		throw(-1);
	}
	assert (AOut->M == (nCells * nVars));
	assert (AOut->N == numModes);

	// transfer rows from regressor to AOut
	for (int i = 0; i < gP.size(); i++) {
		cout << (double)i / gP.size() * 100 << " percent points emplaced \r";
		for (int j = 0; j < nVars; j++)
			AOut->changeContext(regressor, 1, numModes, i + j * gP.size(), 0, gP[i] + j * nCells, 0, false);
	}
	cout << endl;

	if (transA == 'T')
		destroyPMat(regressor, false);

}

