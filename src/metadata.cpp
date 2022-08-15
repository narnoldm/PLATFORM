#include "metadata.hpp"

using namespace ::std;

meta::meta()
{
    isInit = false;
}
meta::~meta()
{
}

meta::meta(vector<string> &iToken)
{

    int offset = 7;

    if (iToken[0] != "input")
    {
        printf("First component of input token should be \"input\", was given as \"%s\"\n", (iToken[0]).c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (iToken[1] != "tecplot")
    {
        printf("Second component of meta input token should be \"binaryset\", was given as \"%s\"\n", (iToken[1]).c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (iToken.size() != offset)
    {
        printf("Meta input token should have %d components, given token has %d components\n", offset, (int)(iToken.size()));
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    token.resize(iToken.size() - offset);
    for (int i = offset; i < iToken.size(); i++)
    {
        token[i - offset] = iToken[i];
    }
    int it0 = stoi(iToken[4]);
    int itf = stoi(iToken[5]);
    int its = stoi(iToken[6]);
    init(it0, itf, its, iToken[2], iToken[3]);
}
meta::meta(int t0, int tf, int ts, string &iPrefix, string &iSuffix)
{
    init(t0, tf, ts, iPrefix, iSuffix);
}

void meta::init(int t0, int tf, int ts, string &iPrefix, string &iSuffix)
{
    snap0 = t0;
    snapF = tf;
    snapSkip = ts;
    prefix = iPrefix;
    suffix = iSuffix;
    if (snapF < snap0)
    {
        printf("snapF (%d) must be greater than or equal to snap0 (%d)\n", snapF, snap0);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    nSets = 0;
    for (int i = snap0; i <= snapF; i = i + snapSkip)
        nSets++;
    cout << prefix + to_string(snap0) + suffix << endl;
    checkSize();
    checkExists();
    isInit = true;
    cout << nPoints << " " << nSets << endl;
}

void meta::checkSize()
{
    FILE *fid;
    string filename = prefix + to_string(snap0) + suffix;
    fid = fopen(filename.c_str(), "rb");
    if (fid == NULL)
    {
        printf("Could not open file at %s\n", filename.c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    int header[2] = {0, 0};
    size_t warn;
    warn = fread(&(header[0]), sizeof(int), 1, fid);
    warn = fread(&(header[1]), sizeof(int), 1, fid);
    assert(header[1] == 1);
    nPoints = header[0];
    fclose(fid);
    return;
}

void meta::checkExists()
{
    FILE *fid;
    int rank;
    string filename;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
    {
        for (int i = snap0; i <= snapF; i = i + snapSkip)
        {
            filename = prefix + to_string(i) + suffix;
            cout << filename << endl;
            fid = fopen(filename.c_str(), "rb");
            if (fid == NULL)
            {
                printf("Could not open file at %s\n", filename.c_str());
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
            fclose(fid);
        }
    }
}

void meta::readSingle(int fileID, double *point)
{
    cout << "meta read " << fileID << endl;
    string filename = prefix + to_string(fileID) + suffix;
    cout << filename << endl;
    FILE *fid;
    fid = fopen((prefix + to_string(fileID) + suffix).c_str(), "rb");
    int header[2] = {0, 0};
    size_t warn;
    warn = fread(&(header[0]), sizeof(int), 1, fid);
    warn = fread(&(header[1]), sizeof(int), 1, fid);
    if (header[0] != nPoints)
    {
        printf("Binary file at %s has leading dimension is %d, expected to be %d\n", filename.c_str(), header[0], (int)nPoints);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (header[1] != 1)
    {
        printf("Binary file at %s has leading dimension is %d, expected to be 1\n", filename.c_str(), header[1]);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    warn = fread(point, sizeof(double), nPoints, fid);
    fclose(fid);
}

void meta::readSingleLowMem(int fileID, pMat* dataMat, int colIdx)
{
    string filename = prefix + to_string(fileID) + suffix;
    readSingleLowMem(filename, dataMat, colIdx);
}

void meta::readSingleLowMem(string filename, pMat* dataMat, int colIdx)
{

    cout << filename << "\r" << flush;
    // cout << "WARNING: meta low-memory read only works if any SZPLT is reordered" << endl;
    // cout << "Make sure activateReorder has been called on any tecIO metadata objects" << endl;

    long maxReadSize = min(nPoints, (long)MAX_IOVEC_SIZE);

    vector<double> getD(maxReadSize);
    int numLoops = ceil((float)nPoints / (float)maxReadSize);

    // open file and read header
    FILE *fid;
    int header[2];
    size_t warn;
    fid = fopen((filename).c_str(), "rb");
    warn = fread(&(header[0]), sizeof(int), 1, fid);
    warn = fread(&(header[1]), sizeof(int), 1, fid);
    if ((header[0] != nPoints) || (header[1] != 1))
    {
        cout << "Unexpected binary shape: " << header[0] << " " << header[1] << endl;
        cout << "Should be: " << nPoints << " 1" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // ingest file
    int readSize, dataIdx;
    for (int i = 0; i < numLoops; ++i)
    {
        if (i < numLoops - 1)
            readSize = maxReadSize;
        else
            readSize = nPoints - maxReadSize * (numLoops - 1);

        warn = fread(&(getD[0]), sizeof(double), readSize, fid);
        for (int j = 0; j < readSize; ++j)
        {
            dataIdx = dataMat->getDataIndex(i * maxReadSize + j, colIdx);
            if (dataIdx >= 0)
            {
                dataMat->dataD[dataIdx] = getD[j];
            }
        }
    }

    fclose(fid);

}

void meta::batchRead(pMat *loadMat)
{
    double t1, t2;
    t1 = MPI_Wtime();
    if (loadMat->mb == nPoints)
    {
        int iP = 0;
        int fileIndex = snap0;
        int localC = 0;
        for (int i = 0; i < nSets; i++)
        {
            iP = (int)(i / loadMat->nb);
            while (iP > (loadMat->pG->pcol - 1))
            {
                iP = iP - loadMat->pG->pcol;
            }
            if (loadMat->pG->rank == iP)
            {
                fileIndex = snap0 + i * snapSkip;
                cout << "Proc " << iP << " is reading file " << fileIndex << endl;
                readSingle(fileIndex, loadMat->dataD.data() + nPoints * localC);
                localC++;
            }
        }
        cout << "Eaiting for other processes read" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        cout << "Even read" << endl;
        for (int j = 0; j < nSets; j++)
        {
            // check if process owns a part of this column
            if (loadMat->pG->mycol == (j / loadMat->nb) % loadMat->pG->pcol)
            {
                readSingleLowMem(snap0 + j * snapSkip, loadMat, j);
            }
        }
    }

    t2 = MPI_Wtime();
    cout << endl << "Batch read took " << t2 - t1 << " secs" << endl;

}

void meta::batchRead(pMat *loadMat, int ii)
{
    if (loadMat->mb == nPoints)
    {
        int iP = 0;
        int fileIndex = snap0;
        int localC = 0;
        int i = ii;
        iP = (int)(i / loadMat->nb);
        while (iP > (loadMat->pG->pcol - 1))
        {
            iP = iP - loadMat->pG->pcol;
        }
        if (loadMat->pG->rank == iP)
        {
            fileIndex = snap0 + i * snapSkip;

            cout << "proc " << iP << " is reading file " << fileIndex << "\r" << flush;

            readSingle(fileIndex, loadMat->dataD.data() + nPoints * localC);
            localC++;
        }
        cout << endl;
        cout << "waiting for other processes read" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        cout << "Even read single" << endl;
        int j = ii;
        if (loadMat->pG->mycol == (j / loadMat->nb) % loadMat->pG->pcol)
        {
            readSingleLowMem(snap0 + j * snapSkip, loadMat, j);
        }
    }
}

void meta::writeSingle(int fileID, double *point, string fpref)
{
    writeSingle(fileID, point, fpref, nPoints, 0);
}

void meta::writeSingle(int fileID, double *point, string fpref, int points, int mode)
{
    // NOTE: mode is irrelevant since this is binary-only

    cout << "meta write single" << endl;
    FILE *fid;
    fid = fopen((fpref + to_string(fileID) + suffix).c_str(), "wb");

    const int ONE = 1;
    fwrite(&points, sizeof(int), 1, fid);
    fwrite(&ONE, sizeof(int), 1, fid);
    fwrite(point, sizeof(double), points, fid);
    fclose(fid);

}

void meta::batchWrite(pMat *loadMat)
{
    batchWrite(loadMat, "out/", prefix);
}
void meta::batchWrite(pMat *loadMat, string dir, string fpref)
{
    batchWrite(loadMat, dir, fpref, 0, nSets, 1);
}
void meta::batchWrite(pMat *loadMat, string dir, string fpref, int nModes)
{
    batchWrite(loadMat, dir, fpref, 0, nModes, 1);
}
void meta::batchWrite(pMat *loadMat, string dir, string fpref, int mStart, int mEnd, int mSkip)
{
    batchWrite(loadMat, dir, fpref, mStart, mEnd, mSkip, snap0, snapSkip, 0, 0);
}

// Write columns or rows of loadMat to individual binary files
// dir: directory to write data to
// fpref: output file prefix
// mStart: starting row/column index of loadMat submatrix to be written (zero-indexed)
// mEnd: ending row/column index of loadMat submatrix to be written (zero-indexed)
// mSkip: row/column index increment of loadMat submatrix to be written
// fStart: starting index of output file names
// fSkip: index increment of output file names
// dim: if 0, write columns of loadMat, otherwise write rows
// mode:
//  if 0, write derived type format and binary
//  if 1, only write derived type format
//  if 2, only write binary
void meta::batchWrite(pMat *loadMat, string dir, string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip, int dim, int mode)
{

    // check if system commands work
    if (!system(NULL))
    {
        printf("Could not access system commands, please diagnose\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (mode > 2)
    {
        cout << "Invalid mode passed to batchWrite: " << mode << endl;
        throw(-1);
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
        int ierr = system(("mkdir " + dir).c_str());

    MPI_Barrier(MPI_COMM_WORLD);
    int iP = 0, fileIndex, localC = 0;

    if (!isInit)
    {
        nPoints = loadMat->M;
        nSets = loadMat->N;
    }

    double t1, t2;
    t1 = MPI_Wtime();
    int vecLength = 0;
    int currentIdx = 0;
    int procRowOrCol, procRowOrColOpp;
    vector<double> tempR;
    MPI_Comm col_comms;
    if (dim == 0)
    {
        vecLength = nPoints;
        loadMat->commCreate(col_comms, 0);
        procRowOrCol = loadMat->pG->mycol;
        procRowOrColOpp = loadMat->pG->myrow;
    }
    else
    {
        vecLength = nSets;
        loadMat->commCreate(col_comms, 1);
        procRowOrCol = loadMat->pG->myrow;
        procRowOrColOpp = loadMat->pG->mycol;
    }
    tempR.resize(vecLength, 0);

    int rowColCheck;
    for (int j = mStart; j < mEnd; j = j + mSkip)
    {
        fill(tempR.begin(), tempR.end(), 0.0);

        // check whether process owns a piece of this column (row)
        if (dim == 0)
        {
            rowColCheck = (j / loadMat->nb) % loadMat->pG->pcol; // process column index
        }
        else
        {
            rowColCheck = (j / loadMat->mb) % loadMat->pG->prow;  // process row index
        }
        if (procRowOrCol == rowColCheck)
        {
            if (dim == 0)
            {
                for (int i = 0; i < nPoints; i++)
                {
                    int xi = i % loadMat->mb;
                    int li = i / (loadMat->pG->prow * loadMat->mb);

                    if (loadMat->pG->myrow == (i / loadMat->mb) % loadMat->pG->prow)
                    {
                        tempR[i] = loadMat->dataD[currentIdx * loadMat->myRC[0] + li * loadMat->mb + xi];
                    }
                }
            }
            else
            {
                for (int i = 0; i < nSets; i++)
                {
                    int xi = currentIdx % loadMat->mb;  // local block row index
                    int li = currentIdx / (loadMat->pG->prow * loadMat->mb); // subblock row index
                    int mi = i / (loadMat->pG->pcol * loadMat->nb); // subblock column index
                    int localcol = mi * loadMat->nb + i % loadMat->nb; // local data column index

                    if (loadMat->pG->mycol == (i / loadMat->nb) % loadMat->pG->pcol)  // process column index
                    {
                        tempR[i] = loadMat->dataD[localcol * loadMat->myRC[0] + li * loadMat->mb + xi];
                    }
                }
            }
            currentIdx++;
        }
        MPI_Allreduce(MPI_IN_PLACE, tempR.data(), tempR.size(), MPI_DOUBLE, MPI_SUM, col_comms);
        // if process owns a piece of vector and is the first process in the process row/column
        if ((procRowOrCol == rowColCheck) && (procRowOrColOpp == 0))
        {
            fileIndex = fStart + (j - mStart) * fSkip;
            printf("proc %d is writing %d\n", loadMat->pG->rank, fileIndex);
            writeSingle(fileIndex, tempR.data(), dir + "/" + fpref, vecLength, mode);
        }
    }
    tempR.clear();
    cout << endl;
    MPI_Comm_free(&col_comms);
    t2 = MPI_Wtime();
    cout << "batch Write took " << t2 - t1 << " secs" << endl;
}

/*
    ===================
    START tecIO METHODS
    ===================
*/

tecIO::tecIO(int t0, int tf, int ts, string &iPrefix, string &iSuffix)
{
    init(t0, tf, ts, iPrefix, iSuffix, "");
}

tecIO::tecIO()
{
    isInit = false;
}

tecIO::tecIO(vector<string> &iToken, string cellIDF)
{

    int offset = 7;
    if (iToken[0] != "input")
    {
        printf("First component of input token should be \"input\", was given as \"%s\"\n", (iToken[0]).c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (iToken[1] != "tecplot")
    {
        printf("Second component of tecIO input token should be \"tecplot\", was given as \"%s\"\n", (iToken[1]).c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    token.resize(iToken.size() - offset);
    for (int i = offset; i < iToken.size(); i++)
    {
        token[i - offset] = iToken[i];
    }
    int it0 = stoi(iToken[4]);
    int itf = stoi(iToken[5]);
    int its = stoi(iToken[6]);
    init(it0, itf, its, iToken[2], iToken[3], cellIDF);
}

void tecIO::init(int t0, int tf, int ts, string &iPrefix, string &iSuffix, string cellIDF)
{
    snap0 = t0;
    snapF = tf;
    snapSkip = ts;
    prefix = iPrefix;
    suffix = iSuffix;
    if (snapF < snap0)
    {
        printf("snapF (%d) must be greater than or equal to snap0 (%d)\n", snapF, snap0);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    nSets = 0;
    varName.clear();
    varIndex.clear();
    for (int i = snap0; i <= snapF; i = i + snapSkip)
        nSets++;
    checkSize();
    isInit = true;

    // set cellIDFile
    if (cellIDF == "")
    {
        // default to first file in sequence if not provided
        cellIDFile = prefix + to_string(snap0) + suffix;
    }
    else
    {
        // TODO: allow to read cellID from binary file?
        if (FILE *file = fopen(cellIDF.c_str(), "r"))
        {
            fclose(file);
            cellIDFile = cellIDF;
        }
        else
        {
            printf("Could not find cellIDF at %s\n", cellIDF.c_str());
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    cout << "tecIO: " << nPoints << " (" << nCells << " * " << numVars << ") x " << nSets << endl;
}

tecIO::~tecIO()
{
}

void tecIO::checkSize()
{
    getDimNodes();
    for (int i = 0; i < token.size(); i = i + 2)
    {
        cout << "token " << token[i] << endl;
        addVar(token[i], token[i + 1]);
    }
    nPoints = nCells * numVars;
}

void tecIO::checkExists()
{
    void *fH = NULL;
    int rank, size;
    string filename;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = snap0; i <= snapF; i = i + snapSkip)
    {
        if (i % size == rank)
        {
            filename = prefix + to_string(i) + suffix;
            tecFileReaderOpen(filename.c_str(), &fH);
            if (fH == NULL)
            {
                printf("Could not open file at %s \n", filename.c_str());
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
            tecFileReaderClose(&fH);
        }
    }
}

// alternative to readSingle that doesn't require a huge vector
void tecIO::readSingleLowMem(int fileID, pMat* dataMat, int colIdx)
{
    string filename = prefix + to_string(fileID) + suffix;
    readSingleLowMem(filename, dataMat, colIdx);
}

void tecIO::readSingleLowMem(string filename, pMat* dataMat, int colIdx)
{

    if (dataMat->M != nPoints)
    {
        cout << "readSingleLowMem only reads to columns for now" << endl;
        throw(-1);
    }

    cout << filename << "\r" << flush;
    void *fH = NULL;
    tecFileReaderOpen(filename.c_str(), &fH);

    int type, dataIdx;
    long maxReadSize = min(nCells, (long)MAX_IOVEC_SIZE);
    vector<float> getF(maxReadSize);
    vector<double> getD(maxReadSize);

    int numCellLoops = ceil((float)nCells / (float)maxReadSize);
    int readSize, cellIdx;

    for (int i = 0; i < numVars; i++)
    {
        tecZoneVarGetType(fH, 1, varIndex[i], &type);
        for (int j = 0; j < numCellLoops; j++)
        {

            if (j < numCellLoops - 1)
                readSize = maxReadSize;
            else
                readSize = nCells - maxReadSize * (numCellLoops - 1);

            if (type == 1)
            {
                tecZoneVarGetFloatValues(fH, 1, varIndex[i], j * maxReadSize + 1, readSize, getF.data());
                copy(getF.begin(), getF.end(), getD.begin());
            }
            else if (type == 2)
            {
                tecZoneVarGetDoubleValues(fH, 1, varIndex[i], j * maxReadSize + 1, readSize, getD.data());
            }

            for (int k = 0; k < readSize; ++k)
            {
                cellIdx = j * maxReadSize + k;
                if (reorder)
                {
                    dataIdx = dataMat->getDataIndex(i * nCells + cellID[cellIdx], colIdx);
                }
                else
                {
                    dataIdx = dataMat->getDataIndex(i * nCells + cellIdx, colIdx);
                }

                if (dataIdx >= 0)
                {
                    dataMat->dataD[dataIdx] = getD[k];
                }
            }
        }
    }

    tecFileReaderClose(&fH);

}

void tecIO::readSingle(int fileID, double *point)
{
    string filename = prefix + to_string(fileID) + suffix;
    readSingle(filename, point);
}

void tecIO::readSingle(string filename, double *point)
{
    cout << filename << "\r" << flush;
    void *fH = NULL;
    tecFileReaderOpen(filename.c_str(), &fH);
    int type;
    if (fH == NULL)
    {
        printf("Could not open file at %s\n", filename.c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    vector<float> get;
    vector<double> temp;

    if (reorder)
    {
        temp.resize(nCells);
    }

    for (int i = 0; i < numVars; i++)
    {
        tecZoneVarGetType(fH, 1, varIndex[i], &type);
        if (type == 1)
        {
            get.resize(nCells);
            tecZoneVarGetFloatValues(fH, 1, varIndex[i], 1, nCells, get.data());
            for (int j = 0; j < nCells; j++)
            {
                point[j + i * nCells] = (double)get[j];
            }
        }
        else if (type == 2)
        {
            tecZoneVarGetDoubleValues(fH, 1, varIndex[i], 1, nCells, &(point[i * nCells]));
        }
        if (reorder)
        {
            // TODO: there has to be a lower memory alternative
            for (int j = 0; j < nCells; j++)
            {
                temp[j] = point[i * nCells + j];
            }
            for (int j = 0; j < nCells; j++)
            {
                point[i * nCells + j] = temp[idx[j]];
            }
        }
    }
    tecFileReaderClose(&fH);
}

void tecIO::batchWrite_bin(pMat *loadMat, string dir, string fpref)
{
    batchWrite_bin(loadMat, dir, fpref, 0, nSets, 1, snap0, snapSkip);
}

// modification of batchWrite that uses idx and MPI_File_write_at_all() to guarantee write order
void tecIO::batchWrite_bin(pMat* dataMat, string dir, string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip)
{
    if ((numVars * nCells) != dataMat->M)
    {
        cout << "numCells * numVars does not match dataMat->M: " << (numVars * nCells) << " " << dataMat->M << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // make output directory
    // check if system commands work
    if (!system(NULL))
    {
        printf("Could not access system commands, please diagnose\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (!rank)
        int ierr = system(("mkdir " + dir).c_str());
    MPI_Barrier(MPI_COMM_WORLD);

    const char* buffer;
    MPI_Offset offset = 0;

    // file info
    MPI_Status status;
    MPI_File fh;
    int fileIndex;
    string filename;

    int bufInt, dataIdx, dofIdx, bufLen;
    double t1 = MPI_Wtime();
    for (int k = mStart; k < mEnd; k = k + mSkip)
    {

        // check if process owns part of this column
        if (dataMat->pG->mycol == (k / dataMat->nb) % dataMat->pG->pcol)
        {
            fileIndex = fStart + (k - mStart) * fSkip;
            filename = dir + "/" + fpref + to_string(fileIndex) + ".bin";
            MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

            // write header
            if (dataMat->pG->myrow == 0)
            {
                bufInt = dataMat->M;
                MPI_File_write_at_all(fh, 0, &bufInt, 1, MPI_INT, &status);
                bufInt = 1;
                MPI_File_write_at_all(fh, 4, &bufInt, 1, MPI_INT, &status);
            }

            cout << "Binary write " << (k + 1) << endl;
            for (int j = 0; j < numVars; ++j)
            {
                for (int i = 0; i < nCells; ++i)
                {
                    if (reorder)
                    {
                        dofIdx = j * nCells + i;
                    }
                    else
                    {
                        dofIdx = j * nCells + idx[i];
                    }

                    // check if process owns part of this row
                    if (dataMat->pG->myrow == (dofIdx / dataMat->mb) % dataMat->pG->prow)
                    {
                        dataIdx = dataMat->getDataIndex(dofIdx, k);

                        // can do a little optimization since the memory is contiguous
                        if (reorder)
                        {
                            // beginning of process block
                            // up to end of block, or end of this variable's data
                            if ((((j * nCells) + i) % dataMat->mb) == 0)
                            {
                                bufLen = min((long)dataMat->mb, nCells - (long)i);
                            }
                            // beginning of variable data
                            // up to end of block
                            else if (i == 0)
                            {
                                bufLen = min((long)dataMat->mb, dataMat->mb - (j * nCells + i) % dataMat->mb);
                            }
                            else
                            {
                                dataIdx = -1;
                            }
                        }
                        else
                        {
                            bufLen = 1;
                        }

                        if (dataIdx >= 0)
                        {
                            offset = 8 + (j * nCells + i) * sizeof(double);
                            MPI_File_write_at_all(fh, offset, &(dataMat->dataD[dataIdx]), bufLen, MPI_DOUBLE, &status);
                        }
                    }
                }
            }
            MPI_File_close(&fh);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "Binary batch write finished in " << MPI_Wtime() - t1 << " seconds" << endl;

}

void tecIO::writeSingle(int fileID, double *point, string fpref)
{
    writeSingle(fileID, point, fpref, nPoints, 0);
}

void tecIO::writeSingle(int fileID, double *point, string fpref, int points, int mode)
{
    if ((mode == 0) || (mode == 1))
    {
        void *infH = NULL;
        void *outfH = NULL;
        string filename = prefix + to_string(fileID) + suffix;
        if (fixedMesh)
        {
            tecFileReaderOpen((meshFile).c_str(), &infH);
        }
        else
        {
            tecFileReaderOpen(filename.c_str(), &infH);
        }
        if (infH == NULL)
        {
            printf("Could not open file at %s\n", filename.c_str());
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        int zoneType;
        long iMax, jMax, kMax;
        tecZoneGetType(infH, 1, &zoneType);
        tecZoneGetIJK(infH, 1, &iMax, &jMax, &kMax);
        if ((zoneType != 5) && (zoneType != 3))
        {
            printf("Zone is weird/not supported\n");
            MPI_Abort(MPI_COMM_WORLD,-1);
        }

        string varstr = "";
        if (dim == 1)
        {
            varstr = "x";
        }
        if (dim == 2)
        {
            varstr = "x,y";
        }
        if (dim == 3)
        {
            varstr = "x,y,z";
        }
        for (int i = 0; i < numVars; i++)
        {
            varstr = varstr + "," + varName[i];
        }

        filename = fpref + to_string(fileID) + suffix;
        tecFileWriterOpen(filename.c_str(), "Code out", varstr.c_str(), 1, 0, 1, NULL, &outfH);
        if (outfH == NULL)
        {
            printf("Could not open file at %s\n", filename.c_str());
            MPI_Abort(MPI_COMM_WORLD, -1);
        }

        vector<int> varTypes(dim + numVars);
        vector<int> valueLoc(dim + numVars);
        vector<int> passive(dim + numVars);
        vector<int> shareVar(dim + numVars);
        for (int i = 0; i < (dim + numVars); i++)
        {
            if (i < dim)
            {
                tecZoneVarGetType(infH, 1, i + 1, varTypes.data() + i);
                tecZoneVarGetValueLocation(infH, 1, i + 1, valueLoc.data() + i);
                tecZoneVarIsPassive(infH, 1, i + 1, passive.data() + i);
                tecZoneVarGetSharedZone(infH, 1, i + 1, shareVar.data() + i);
            }
            else
            {
                varTypes[i] = 2;
                valueLoc[i] = 0;
                passive[i] = 0;
                shareVar[i] = 0;
            }
        }
        int shareCon, fNeigh, outZone;
        tecZoneConnectivityGetSharedZone(infH, 1, &shareCon);
        tecZoneFaceNbrGetMode(infH, 1, &fNeigh);
        tecZoneCreateFE(outfH, to_string(fileID).c_str(), zoneType, iMax, jMax, &varTypes[0], &shareVar[0], &valueLoc[0], &passive[0], shareCon, 0, 0, &outZone);
        tecZoneSetUnsteadyOptions(outfH, outZone, fileID, (int)fileID);

        // copy nodal coordinate data from baseline file
        vector<float> nfDat;
        vector<double> ndDat;
        for (int i = 0; i < dim; i++)
        {
            if (varTypes[i] == 1)
            {
                nfDat.resize(iMax);
                tecZoneVarGetFloatValues(infH, 1, i + 1, 1, iMax, nfDat.data());
                tecZoneVarWriteFloatValues(outfH, 1, i + 1, 0, iMax, nfDat.data());
            }
            else if (varTypes[i] == 2)
            {
                ndDat.resize(iMax);
                tecZoneVarGetDoubleValues(infH, 1, i + 1, 1, iMax, &ndDat[0]);
                tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, iMax, &ndDat[0]);
            }
        }
        vector<float>().swap(nfDat);
        vector<double>().swap(ndDat);

        // write field data
        vector<double> temp(jMax, 0.0);
        for (int i = dim; i < (dim + numVars); i++)
        {
            if (reorder)
            {
                for (int n = 0; n < jMax; n++)
                {
                    temp[idx[n]] = point[(i - dim) * jMax + n];
                }
                tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, jMax, temp.data());
            }
            else
            {
                tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, jMax, &point[(i - dim) * jMax]);
            }
        }
        vector<double>().swap(temp);

        // write nodemap
        long numValues;
        tecZoneNodeMapGetNumValues(infH, 1, jMax, &numValues);
        vector<int> nodeMap(numValues);
        tecZoneNodeMapGet(infH, 1, 1, jMax, nodeMap.data());
        tecZoneNodeMapWrite32(outfH, 1, 0, 1, numValues, nodeMap.data());
        vector<int>().swap(nodeMap);

        // close files
        tecFileReaderClose(&infH);
        tecFileWriterClose(&outfH);
    }

    if ((mode == 0) || (mode == 2))
    {
        cout << "bin write single " << fileID << endl;
        FILE *fid;
        fid = fopen((fpref + to_string(fileID) + ".bin").c_str(), "wb");

        int ONE = 1;
        fwrite(&nPoints, sizeof(int), 1, fid);
        fwrite(&ONE, sizeof(int), 1, fid);
        for (int i = 0; i < numVars; i++)
        {
            if (reorder)
            {
                for (int j = 0; j < nCells; j++)
                    fwrite(&point[i * nCells + j], sizeof(double), 1, fid);
            }
            else
            {
                for (int j = 0; j < nCells; j++)
                    fwrite(&point[i * nCells + idx[j]], sizeof(double), 1, fid);
            }
        }
        fclose(fid);
    }
}

void tecIO::write_ascii(pMat* dataMat, string filename, string header)
{
    // TODO: generalize
    if (dataMat->N > 1)
    {
        printf("write_ascii only works with column pMats right now\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if (dataMat->block != 0)
    {
        printf("write_ascii only works with square block pMats\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    const char* buffer;
    MPI_Offset offset = 0;

    // open file
    MPI_Status status;
    MPI_File fh;
    MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);

    // write header
    header += "\n";
    if (rank == 0)
    {
        buffer = header.c_str();
        MPI_File_write_at_all(fh, offset, buffer, strlen(buffer), MPI_CHAR, &status);
    }

    // determine length of single line
    ostringstream outstream;
    // account for negative values, possible 3rd exponent digit, add padding where necessary
    double val = -1.0e-200;
    outstream << scientific << setprecision(numeric_limits<double>::digits10) << val << "\n";
    int valLength = strlen((outstream.str()).c_str());
    outstream.str("");
    outstream.clear();

    // TODO: this is all meaningless for general matrix write
    int dataIdx, padding, bufLen;
    double absval;
    if (dataMat->pG->mycol == 0)
    {
        int xi, li;
        string tail;
        string outstring;

        for (int j = 0; j < numVars; ++j)
        {
            for (int i = 0; i < nCells; ++i)
            {

                // can do a little optimization since the memory is contiguous
                if (reorder)
                {
                    // beginning of process block
                    // up to end of block, or end of this variable's data
                    if ((((j * nCells) + i) % dataMat->mb) == 0)
                    {
                        dataIdx = dataMat->getDataIndex(j * nCells + i, 0);
                        bufLen = min((long)dataMat->mb, nCells - (long)i);
                    }
                    // beginning of variable data
                    // up to end of block
                    else if (i == 0)
                    {
                        dataIdx = dataMat->getDataIndex(j * nCells + i, 0);
                        bufLen = min((long)dataMat->mb, dataMat->mb - (j * nCells + i) % dataMat->mb);
                    }
                    else
                    {
                        dataIdx = -1;
                    }
                }
                else
                {
                    dataIdx = dataMat->getDataIndex(j * nCells + idx[i], 0);
                    bufLen = 1;
                }

                if (dataIdx >= 0)
                {
                    padding = 0;
                    tail = "";

                    // add data to stream
                    for (int k = 0; k < bufLen; ++k)
                    {
                        val = dataMat->dataD[dataIdx + k];

                        // pad with additional space if necessary
                        if (val > 1e-100)
                        {
                            padding = 2;
                        }
                        else
                        {
                            if (val == 0.0)
                            {
                                padding = 2;
                            }
                            else
                            {
                                if (val < -1e-100)
                                    padding = 1;
                                else
                                    padding = 0;
                            }
                        }
                        for (int l = 0; l < padding; ++l)
                            tail += " ";
                        tail += "\n";

                        outstream << scientific << setprecision(numeric_limits<double>::digits10) << val << tail;

                    }

                    // write to file
                    outstring = outstream.str();
                    offset = strlen(header.c_str()) + (j * nCells + i) * valLength;
                    buffer = outstring.c_str();
                    MPI_File_write_at_all(fh, offset, buffer, strlen(buffer), MPI_CHAR, &status);
                    outstream.str("");
                    outstream.clear();
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File_close(&fh);

}

void tecIO::addVar(string var, string &norm)
{
    varName.push_back(var);
    int varI = getVariableIndex(var, prefix + to_string(snap0) + suffix);
    varIndex.push_back(varI);
    normID.push_back(norm);
    double temp;
    temp = stod(norm);
    scalingInput.push_back(temp);
    if (varName.size() != varIndex.size())
    {
        printf("Something went wrong with registering variable \"%s\". Was varIndex or varName modified inappropriately?\n", var.c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    numVars = varName.size();
    cout << "metadata has registered " << numVars << " variable(s)" << endl;
}

int tecIO::getVariableIndex(string var, string file)
{
    void *fH;
    int fileVars;
    int tecIndex = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0)
    {

        char *vName = NULL;
        tecFileReaderOpen(file.c_str(), &fH);
        if (fH == NULL)
        {
            printf("Could not open file at %s\n", file.c_str());
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        tecDataSetGetNumVars(fH, &fileVars);
        for (int i = 1; i <= fileVars; i++)
        {
            tecVarGetName(fH, i, &vName);
            if (var == vName)
            {
                printf("%s found in index %d\n", vName, i);
                tecIndex = i;
                break;
            }
            delete[] vName;
            vName = NULL;
        }
        delete[] vName;
        vName = NULL;
        tecFileReaderClose(&fH);
        if (tecIndex == 0)
        {
            cout << "Var not found :" << var << endl;
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
    }
    MPI_Bcast(&tecIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return (tecIndex);
}

void tecIO::getDimNodes()
{

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
    {
        if (fixedMesh)
        {
            checkMeshDim(meshFile);
        }
        else
        {
            checkMeshDim((prefix + to_string(snap0) + suffix));
        }
    }
    MPI_Bcast(&nCells, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
}
void tecIO::checkMeshDim(string filename)
{
    void *fH = NULL;
    long iMax = 0, jMax = 0, kMax = 0;
    long check;
    int id = 0;
    cout << filename << endl;
    tecFileReaderOpen(filename.c_str(), &fH);
    if (fH == NULL)
    {
        printf("Could not open file at %s\n", filename.c_str());
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    tecZoneGetIJK(fH, 1, &iMax, &jMax, &kMax);
    cout << "iMax:" << iMax << ", jMax " << jMax << ", kMax " << kMax << endl;
    check = iMax;
    while (check == iMax)
    {
        id++;
        tecZoneVarGetNumValues(fH, 1, id, &check);
        cout << "var " << id << " has " << check << " values" << endl;
    }
    nCells = jMax;
    printf("var %d is the first cellcentered variable\n", id);
    dim = id - 1;
    printf("dimension is %d\n", dim);
    tecFileReaderClose(&fH);
}

void tecIO::genHash()
{
    genHash("cellIDFile");
}

void tecIO::genHash(string filename)
{
    if (idx.size() != 0)
        return;

    idx.resize(nCells, 0);
    iota(idx.begin(), idx.end(), 0);
    cellID.resize(nCells);
    void *fH;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // default to tecIO.cellIDFile if no argument provided to genHash
    if (filename == "cellIDFile")
    {
        filename = cellIDFile;
    }

    if (filename != "")
    {
        int var_index = 0;
        var_index = getVariableIndex("cell_id", filename);

        // faster if every process loads cell_id instead of broadcasting
        int hashType;
        tecFileReaderOpen(filename.c_str(), &fH);
        if (fH == NULL)
        {
            printf("Could not open file at %s\n", filename.c_str());
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        tecZoneVarGetType(fH, 1, var_index, &hashType);
        if (hashType == 3)
        {
            tecZoneVarGetInt32Values(fH, 1, var_index, 1, nCells, cellID.data());
        }
        else
        {
            printf("genHash() can only read int32 cell_id values\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        tecFileReaderClose(&fH);
        for (int i = 0; i < nCells; i++)
        {
            cellID[i]--;
        }

        stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return cellID[i] < cellID[j]; });
        cout << "Hash table built" << endl;
    }
    else
    {
        cout << "Assume default hash" << endl;
        for (int i = 0; i < nCells; i++)
        {
            cellID[i] = i;
            idx[i] = i;
        }
    }
}

// compute centering profile
// centerMethod = "avg" : centers about the dataset average
// centerMethod = pathToFile : centers about the profile provided by pathToFile
void tecIO::calcCentering(pMat *dataMat, string centerMethod)
{
    calcCentering(dataMat, centerMethod, false);
}

void tecIO::calcCentering(pMat *dataMat, string centerMethod, bool isField)
{
    calcCentering(dataMat, centerMethod, isField, true);
}

void tecIO::calcCentering(pMat *dataMat, string centerMethod, bool isField, bool writeToDisk)
{

    isCentered = true;
    genHash();
    centerVec = new pMat(nPoints, 1, dataMat->pG, 0, 0, 0.0, false);

    cout << "Loading/calculating centering" << endl;

    // update this check when a new method is added
    if ((centerMethod != "avg") && (centerMethod != "avgmag"))
    {

        // load from file
        // TODO: allow reading from small vector
        isField = true;
        if (centerMethod.substr(centerMethod.size()-6, 6) == ".szplt")
        {
            readSingleLowMem(centerMethod, centerVec, 0);
            cout << endl;
        }
        else
        {
            printf("Invalid centering file: %s\n", centerMethod.c_str());
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        printf("Centering calculation has not been re-implemented\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // update global param
    centeringIsField = isField;

    // expand constants to full field for convenience sake
    if (!isField)
    {
        printf("Distributing scalar centering values has not been re-implemented\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // write centering field to file
    cout << "Centering calculated" << endl;
    if (writeToDisk)
    {
        // TODO: get SZPLT to output correctly, without silly fileID requirement
        write_ascii(centerVec, "centerProf.dat", "centerProf");
        cout << "Centering files written" << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

}

void tecIO::centerData(pMat *dataMat)
{
    centerData(dataMat, false);
}

void tecIO::centerData(pMat *dataMat, bool uncenter)
{
    if (!centerVec)
    {
        printf("Centering isn't setup, call calcCentering first\n");
        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    double alpha;
    if (uncenter)
    {
        cout << "Uncentering" << endl;
        alpha = 1.0;
    }
    else
    {
        cout << "Centering" << endl;
        alpha = -1.0;
    }

    // subtract centering vector, column-by-column
    for (int j = 0; j < dataMat->N; ++j)
    {
        dataMat->matrix_Sum('N', dataMat->M, 1, centerVec, 0, 0, 0, j, alpha, 1.0);
    }

    if (uncenter)
    {
        cout << "Data uncentered" << endl;
    }
    else
    {
        cout << "Data centered" << endl;
    }

}

void tecIO::calcScaling(pMat* dataMat, string scaleMethod)
{
    calcScaling(dataMat, scaleMethod, false);
}

void tecIO::calcScaling(pMat* dataMat, string scaleMethod, bool isField)
{
    calcScaling(dataMat, scaleMethod, isField, true);
}

void tecIO::calcScaling(pMat *dataMat, string scaleMethod, bool isField, bool writeToDisk)
{

    isScaled = true;
    genHash();
    scalingSubVec = new pMat(nPoints, 1, dataMat->pG, 0, 0, 0.0, false);
    scalingDivVec = new pMat(nPoints, 1, dataMat->pG, 0, 0, 0.0, false);
    scalingSubVecFull = new pMat(nPoints, 1, dataMat->pG, 0, 0, 0.0, false);

    if ((scaleMethod == "sphere") && (isField))
    {
        printf("If choosing scaleMethod = sphere, you must set scaleIsField = false\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    if ((scaleMethod == "standardize") && (!isField))
    {
        cout << "WARNING: scalar standardization doesn't really make sense." << endl;
    }

    // check if the user has specified ANY scaling factors
    // if so, all scaling will be by scalar values
    for (int i = 0; i < numVars; ++i)
    {
        if (scalingInput[i] > 0.0)
        {
            if (!isField)
            {
                cout << "WARNING: Passed a scaling factor, all scaling will be by scalars" << endl;
            }
            isField = false;
            break;
        }
    }

    // check if any values need to be calculated
    bool calc = false;
    for (int i = 0; i < numVars; ++i)
    {
        if (scalingInput[i] < 0.0)
        {
            calc = true;
            break;
        }
    }

    // update this check when a new method is added
    if ((scaleMethod != "minmax") &&
        (scaleMethod != "standardize") &&
        (scaleMethod != "sphere"))
    {
        // load from file
        // TODO: allow reading from ASCII file
        isField = true;
        ifstream szlSub((scaleMethod + "SubProf.szplt").c_str());
        ifstream szlDiv((scaleMethod + "DivProf.szplt").c_str());
        if (szlSub.good() && szlDiv.good())
        {
            readSingleLowMem(scaleMethod + "SubProf.szplt", scalingSubVec, 0);
            cout << endl;
            readSingleLowMem(scaleMethod + "DivProf.szplt", scalingSubVec, 0);
            cout << endl;
        }
        else
        {
            printf("Could not find scaling files with prefix \"%s\"\n",scaleMethod.c_str());
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        if (calc)
        {
            printf("Scaling calculation has not been re-implemented\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }

    // update global param
    scalingIsField = isField;

    // expand constants to full field
    if (!isField)
    {
        for (int j = 0; j < numVars; ++j)
        {
            for (int i = 0; i < nCells; ++i)
            {
                scalingDivVec->setElement(j * nCells + i, 0, scalingInput[j], false);
            }
        }
    }

    // prevent divisive factors close to zero
    // really only happens when min and max are same, implying field is uniform
    for (int i = 0; i < scalingDivVec->dataD.size(); ++i)
    {
        if (abs(scalingDivVec->dataD[i]) < 1e-15)
        {
            scalingDivVec->dataD[i] = 1.0;
        }
    }

    // write scaling fields to file
    cout << "Scaling calculated" << endl;
    if (writeToDisk)
    {
        // TODO: get SZPLT to output correctly, without silly fileID requirement
        write_ascii(scalingSubVec, "scalingSubProf.dat", "scalingSubProf");
        write_ascii(scalingDivVec, "scalingDivProf.dat", "scalingDivProf");

        // if also centering, combine into single profile (useful for GEMS)
        if (isCentered)
        {
            scalingSubVecFull->matrix_Sum('N', nPoints, 1, centerVec, 0, 0, 0, 0, 1.0, 1.0);
            scalingSubVecFull->matrix_Sum('N', nPoints, 1, scalingSubVec, 0, 0, 0, 0, 1.0, 1.0);
            cout << "Writing combined subtractive scaling factors" << endl;
            write_ascii(scalingSubVecFull, "scalingSubProfFull.dat", "scalingSubProfFull");
        }

        cout << "Scaling files written" << endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

void tecIO::scaleData(pMat *dataMat)
{
    scaleData(dataMat, false);
}

void tecIO::scaleData(pMat *dataMat, bool unscale)
{

    if (!scalingDivVec || !scalingSubVec)
    {
        printf("Scaling isn't setup, call calcScaling first\n");
        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    if (unscale)
    {
        cout << "Unscaling" << endl;
    }
    else
    {
        cout << "Scaling" << endl;
    }

    double divFac;
    if (unscale)
    {
        // multiply rows
        for (int j = 0; j < dataMat->N; ++j)
        {
            dataMat->matrix_elem_mult('N', dataMat->M, 1, 1.0, scalingDivVec, 0, 0, 0, j);
        }

        // add columns
        for (int j = 0; j < dataMat->N; ++j)
        {
            dataMat->matrix_Sum('N', dataMat->M, 1, scalingSubVec, 0, 0, 0, j, 1.0, 1.0);
        }
    }
    else
    {
        // subtract columns
        for (int j = 0; j < dataMat->N; ++j)
        {
            dataMat->matrix_Sum('N', dataMat->M, 1, scalingSubVec, 0, 0, 0, j, -1.0, 1.0);
        }

        // change scalingDivVec to multiplicative factor
        for (int j = 0; j < scalingDivVec->dataD.size(); ++j)
        {
            scalingDivVec->dataD[j] = 1.0 / scalingDivVec->dataD[j];
        }

        // divide rows
        for (int j = 0; j < dataMat->N; ++j)
        {
            dataMat->matrix_elem_mult('N', dataMat->M, 1, 1.0, scalingDivVec, 0, 0, 0, j);
        }

        // reverse scalingDivVec modification
        for (int j = 0; j < scalingDivVec->dataD.size(); ++j)
        {
            scalingDivVec->dataD[j] = 1.0 / scalingDivVec->dataD[j];
        }

    }

    if (unscale)
    {
        cout << "Data unscaled" << endl;
    }
    else
    {
        cout << "Data scaled" << endl;
    }

}

// compute cell-wise quantity across all snapshots
// methodName: accepts "avg", "avgmag", "min", "max", "l2"
void tecIO::calcGroupQuant(pMat *dataMat, double &outVal, vector<double> &outVec, int varIdx, string methodName, bool outField)
{
    double val, groupVal;
    vector<double> valVec(nCells);

    // prefill aggregate vector
    if (methodName == "min") {
        fill(valVec.begin(), valVec.end(), HUGE_VAL);
    }
    else if (methodName == "max")
    {
        fill(valVec.begin(), valVec.end(), -HUGE_VAL);
    }
    else
    {
        fill(valVec.begin(), valVec.end(), 0.0);
        if (methodName == "l2")
        {
            outVal = -HUGE_VAL;
        }
    }

    // loop snapshots
    for (int j = 0; j < nSets; ++j)
    {

        if (methodName == "l2")
        {
            fill(outVec.begin(), outVec.end(), 0.0);
        }

        for (int i = 0; i < nCells; ++i)
        {
            if (methodName == "min")
            {
                groupVal = HUGE_VAL;
            }
            else if (methodName == "max")
            {
                groupVal = -HUGE_VAL;
            }
            else
            {
                groupVal = 0.0;
            }

            // get contributions from all fields in group
            for (int g = 0; g < numVars; ++g)
            {
                if (scalingInput[g] == scalingInput[varIdx])
                {
                    // gather scaling values
                    // expand as necessary for new methods
                    if (methodName == "min")
                    {
                        val = dataMat->getLocalElement(g * nCells + i, j, HUGE_VAL);
                        if (val != HUGE_VAL)
                        {
                            groupVal = min(val, groupVal);
                        }
                    }
                    else if (methodName == "max")
                    {
                        val = dataMat->getLocalElement(g * nCells + i, j, -HUGE_VAL);
                        if (val != -HUGE_VAL)
                        {
                            groupVal = max(val, groupVal);
                        }
                    }
                    else if (methodName == "avg")
                    {
                        val = dataMat->getLocalElement(g * nCells + i, j);
                        groupVal += val;
                    }
                    else if (methodName == "avgmag")
                    {
                        // MUST use global getElement, proc may not own all variables in group
                        val = dataMat->getElement(g * nCells + i, j);
                        groupVal += pow(val, 2);
                    }
                    else if (methodName == "l2")
                    {
                        val = dataMat->getLocalElement(g * nCells + i, j);
                        groupVal += pow(val, 2);
                    }
                }
            }
            if (methodName == "min")
            {
                valVec[i] = min(valVec[i], groupVal);
            }
            else if (methodName == "max")
            {
                valVec[i] = max(valVec[i], groupVal);
            }
            else if (methodName == "avg")
            {
                valVec[i] += groupVal;
            }
            else if (methodName == "avgmag")
            {
                valVec[i] += sqrt(groupVal);
            }
            else if (methodName == "l2")
            {
                valVec[i] = groupVal;
            }
        }

        // compute l2 norm for this snapshot, take max
        if (methodName == "l2")
        {
            MPI_Allreduce(valVec.data(), outVec.data(), nCells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
            val = 0.0;
            for (int i = 0; i < nCells; ++i)
            {
                val += outVec[i];
            }
            val = sqrt(val);
            outVal = max(outVal, val);
        }

    }

    // only collected local elements, so need to reduce
    if (methodName == "min")
    {
        MPI_Allreduce(valVec.data(), outVec.data(), nCells, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    }
    else if (methodName == "max")
    {
        MPI_Allreduce(valVec.data(), outVec.data(), nCells, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    }
    else if (methodName == "avg")
    {
        MPI_Allreduce(valVec.data(), outVec.data(), nCells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int i = 0; i < nCells; ++i)
        {
            outVec[i] /= nSets;
        }
    }
    else if (methodName == "avgmag")
    {
        // don't need to reduce because this used global getElement
        for (int i = 0; i < nCells; ++i)
        {
            outVec[i] = valVec[i] / nSets;
        }
    }
    else if (methodName == "l2")
    {
        // already reduced
        for (int i = 0; i < nCells; ++i)
        {
            outVec[i] = outVal;
        }
    }

    if (!outField)
    {
        if (methodName == "min")
        {
            outVal = HUGE_VAL;
            for (int i = 0; i < nCells; ++i)
            {
                outVal = min(outVal, outVec[i]);
            }
        }
        else if (methodName == "max")
        {
            outVal = -HUGE_VAL;
            for (int i = 0; i < nCells; ++i)
            {
                outVal = max(outVal, outVec[i]);
            }
        }
        else if ((methodName == "avg") || (methodName == "avgmag"))
        {
            outVal = 0.0;
            for (int i = 0; i < nCells; ++i)
            {
                outVal += outVec[i];
            }
            outVal /= nCells;
        }
    }

}

void tecIO::activateReorder()
{
    activateReorder("cellIDFile");
}

void tecIO::activateReorder(string file)
{
    genHash(file);
    reorder = true;
}

// check if metadata for two meta objects are the same, EXCEPT for file numbers
int compareMeta(meta* meta1, meta* meta2)
{

	// check path members
	if ( (meta1->prefix == meta2->prefix) && (meta1->suffix == meta2->suffix) )
    {

		if ( (meta1->snap0 == meta2->snap0) &&
			 (meta1->snapF == meta2->snapF) &&
			 (meta1->snapSkip == meta2->snapSkip) )
        {
			// objects are identical
			return(1);
		}
        else
        {
			// objects differ only in the file counts
			return(2);
		}
		return(true);

	} else {
		// objects are totally different
		return(0);
	}

}
