#include "sampling.hpp"

using namespace ::std;

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

void eigenvector_oversampling(pMat* URes, pMat* USol, int sampMethod, int nCells, int nVars, int nDOF, int numModesRHS,
							  int PointsNeeded, set<int>& samplingPoints, vector<int>& gP, string& timingOutput) {

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// declare some variables for this routine
	pMat *U_E, *VT_E;
	vector<double> S_E;
	int M, N, minDim;
	int numInitPoints, numCurrentPoints, numCurrentDOFs;
	int cellID;

	PGrid* evenG = URes->pG; // reuse process grid

	// greedy sampling metric vector
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
	pMat *URHS_samp_E = new pMat(PointsNeeded*nVars, numModesRHS, evenG, false);
	pMat *URHS_samp_E_copy = new pMat(PointsNeeded*nVars, numModesRHS, evenG, false);
	for (int j = 0; j < numCurrentPoints; ++j) {
		cout << (double)j / numCurrentPoints * 100 << " percent points extracted \r";
		for (int k = 0; k < nVars; ++k)
			URHS_samp_E_copy->changeContext(URes, 1, numModesRHS, gP[j] + k * nCells, 0, j * nVars + k, 0, false);
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
		rVec->matrix_Product('N', 'T', 1, nDOF, minDim, VT_E, minDim-1, 0, URes, 0, 0, 1.0, 0.0, 0, 0);
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
			URHS_samp_E_copy->changeContext(URes, 1, numModesRHS, k * nCells + gP.back(), 0, (numCurrentPoints-1) * nVars + k, 0, false);
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
}

void classic_greedy_oversampling(pMat* URes, pMat* USol, int sampMethod, int nCells, int nVars, int nDOF, int numModesRHS,
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

	// variables for computing argsort
	int breakSignal;
	double maxLocal, maxGlobal;
	int argMaxLocal, argMaxGlobal;

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
}

