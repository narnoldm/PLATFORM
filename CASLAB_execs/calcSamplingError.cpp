// Compute the l2-norm of the pseudo-inverse of a subsampled basis
// This is equivalent to the maximum singular value of the pseudo-inverse of the subsampled basis

#include "metadata.hpp"
#include "param.hpp"

using namespace ::std;

int main(int argc, char *argv[])
{
	// some setup
	MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();
	PGrid *evenG;
    evenG = new PGrid(rank, size, 0);
	double t1, t2;

	// get input file
	paramMap inputFile("sampErr.inp", rank);

	// specify rank to write to STDOUT
	int debug_proc = 0;
	inputFile.getParamInt("stdout_proc", debug_proc);
    if (rank != debug_proc)
    {
        cout.rdbuf(sink.rdbuf());
    }

	// prep basis dataset
	string inputBasis;
	inputFile.getParamString("basisInputString", inputBasis); 	// PDP-demarcated input for basis
	cout << "Basis input string is: " << inputBasis << endl;
	vector<string> tokenBasis;
	tokenparse(inputBasis, "|", tokenBasis);
	meta *datasetBasis = new meta(tokenBasis);
	int nModes = datasetBasis->nSets;

	// have to provide number of cells and variables, since reading from binary (not SZPLT)
	int nCells, nVars;
	inputFile.getParamInt("nCells", nCells); 
	inputFile.getParamInt("nVars", nVars);

	// directory where samplingPoints.txt is stored, and output file will be written
	string pinvDir;
	inputFile.getParamString("pinvDir", pinvDir);

	// read modes from disk
	cout << "Loading basis..." << endl;
	pMat *UIn = new pMat(datasetBasis->nPoints, datasetBasis->nSets, evenG, false);
	datasetBasis->batchRead(UIn);
	
	// load sampling file
	cout << "Getting sample indices..." << endl;
	string sampleFileName = pinvDir + "/samplingPoints.txt";
	fstream sampFile(sampleFileName, ios_base::in);
	int nSamps;
	sampFile >> nSamps;
	vector<int> sampleIdxs(nSamps);
	for (int i = 0; i < nSamps; ++i) {
		sampFile >> sampleIdxs[i];
		--sampleIdxs[i]; 
	}
	sampFile.close();

	// extract rows of basis
	cout << "Extracting basis rows..." << endl;
	t1 = MPI_Wtime();
	pMat *UIn_samp = new pMat(nSamps * nVars, nModes, evenG, false);
    for (int i = 0; i < nSamps; i++) {
        cout << (double)i / nSamps * 100 << " percent points extracted \r";
		for (int j = 0; j < nVars; j++)
			UIn_samp->changeContext(UIn, 1, nModes, sampleIdxs[i] + j * nCells, 0, i + j * nSamps, 0, false);
    }
    t2 = MPI_Wtime();
    cout << endl << "Extraction of basis rows took " << t2 - t1 << " seconds" << endl;
	destroyPMat(UIn, false);

	// compute pseudo-inverse of sampled basis
	cout << "Computing pseudo-inverse..." << endl;
    pMat *pinvUIn_samp = new pMat(UIn_samp->N, UIn_samp->M, evenG, false);
    pinvUIn_samp->pinv(UIn_samp);
	destroyPMat(UIn_samp, false);
	
	// compute SVD of pseudo-invers of sampled basis
	cout << "Computing SVD..." << endl;
	int minDim = min(pinvUIn_samp->M, pinvUIn_samp->N);
	vector<double> S(minDim);
	pMat *U = new pMat(pinvUIn_samp->M, minDim, evenG, false);
	pMat *VT = new pMat(minDim, pinvUIn_samp->N, evenG, false);

	pinvUIn_samp->svd_run(pinvUIn_samp->M, pinvUIn_samp->N, 0, 0, U, VT, S, false);
	destroyPMat(pinvUIn_samp, false);
	destroyPMat(U, false);
	destroyPMat(VT, false);

	// get maximum singular value, this is the error measure
	if (rank == 0) {
		cout << "Maximum singular value of [P^T U]^+: " << setprecision(16) << S[0] << endl;
		string outputFileName = pinvDir + "/samplingError.dat";
		ofstream outputFile(outputFileName);
		outputFile << setprecision(16) << S[0] << endl;
		outputFile.close();
	}

	// finalize
	cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}