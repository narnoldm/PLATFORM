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

    assert(iToken[0] == "input");
    assert(iToken[1] == "binaryset");
    assert(iToken.size() == offset);
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
    assert(snapF >= snap0);
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
    fid = fopen((prefix + to_string(snap0) + suffix).c_str(), "rb");
    assert(fid != NULL);
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
    {
        for (int i = snap0; i <= snapF; i = i + snapSkip)
        {
            cout << (prefix + to_string(i) + suffix) << endl;
            fid = fopen((prefix + to_string(i) + suffix).c_str(), "rb");
            assert(fid != NULL);
            fclose(fid);
        }
    }
}

void meta::readSingle(int fileID, double *point)
{
    cout << "meta read " << fileID << endl;
    cout<< (prefix + to_string(fileID) + suffix)<<endl;
    FILE *fid;
    fid = fopen((prefix + to_string(fileID) + suffix).c_str(), "rb");
    int header[2] = {0, 0};
    size_t warn;
    warn = fread(&(header[0]), sizeof(int), 1, fid);
    warn = fread(&(header[1]), sizeof(int), 1, fid);
    assert(header[1] == 1);
    assert(header[0] == nPoints);
    warn = fread(point, sizeof(double), nPoints, fid);
    fclose(fid);
}

void meta::readSingleLowMem(int fileID, pMat* dataMat, int colIdx)
{
    cout << "No low-memory call for meta" << endl;
    throw(-1);
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
    writeSingle(fileID, point, fpref, nPoints);
}

void meta::writeSingle(int fileID, double *point, string fpref, int points)
{

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
    batchWrite(loadMat, dir, fpref, mStart, mEnd, mSkip, snap0, snapSkip, 0);
}

// Write columns or rows of loadMat to individual binary files
// dir: directory to write data to
// fpref: output file prefix
// mStart: starting row/column index of loadMat submatrix to be written (zero-indexed)
// mEnd: ending row/column index of loadMat submatrix to be written (zero-indexed)
// mSkip: row/column index increment of loadMat submatrix to be written
// fStart: starting index of output file names
// fSkip: index increment of output file names
// writeCols: if true, write columns of loadMat, otherwise write rows
void meta::batchWrite(pMat *loadMat, string dir, string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip, int dim)
{

    assert(system(NULL)); //check if system commands work
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
        // cout << j << " " << procRowOrCol << " " << rowColCheck << " " << procRowOrColOpp << endl;
        if ((procRowOrCol == rowColCheck) && (procRowOrColOpp == 0))
        {
            fileIndex = fStart + (j - mStart) * fSkip;
            printf("proc %d is writing %d\n", loadMat->pG->rank, fileIndex);
            writeSingle(fileIndex, tempR.data(), dir + "/" + fpref, vecLength);
        }
    }
    tempR.clear();
    cout << endl;
    MPI_Comm_free(&col_comms);
    t2 = MPI_Wtime();
    cout << "batch Write took " << t2 - t1 << " secs" << endl;
}

tecIO::tecIO(int t0, int tf, int ts, string &iPrefix, string &iSuffix)
{
    init(t0, tf, ts, iPrefix, iSuffix);
}
tecIO::tecIO()
{
    isInit = false;
}
tecIO::tecIO(vector<string> &iToken)
{

    int offset = 7;
    assert(iToken[0] == "input");
    assert(iToken[1] == "tecplot");

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
void tecIO::init(int t0, int tf, int ts, string &iPrefix, string &iSuffix)
{
    snap0 = t0;
    snapF = tf;
    snapSkip = ts;
    prefix = iPrefix;
    suffix = iSuffix;
    assert(snapF >= snap0);
    nSets = 0;
    varName.clear();
    varIndex.clear();
    for (int i = snap0; i <= snapF; i = i + snapSkip)
        nSets++;
    checkSize();
    //checkExists();
    isInit = true;
    cout << nPoints << " " << nSets << endl;
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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    for (int i = snap0; i <= snapF; i = i + snapSkip)
    {
        if (i % size == rank)
        {
            tecFileReaderOpen((prefix + to_string(i) + suffix).c_str(), &fH);
            assert(fH != NULL);
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
                    dataIdx = dataMat->getDataIndex(i * nCells + idx[cellIdx], colIdx);
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
    assert(fH != NULL);
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

void tecIO::writeSingle(int fileID, double *point, string fpref)
{
    writeSingle(fileID, point, fpref, nPoints);
}

void tecIO::writeSingle(int fileID, double *point, string fpref, int points)
{
    void *infH = NULL;
    void *outfH = NULL;
    if (fixedMesh)
    {
        tecFileReaderOpen((meshFile).c_str(), &infH);
    }
    else
    {
        tecFileReaderOpen((prefix + to_string(fileID) + suffix).c_str(), &infH);
    }
    assert(infH != NULL);
    int zoneType;
    long iMax, jMax, kMax;
    tecZoneGetType(infH, 1, &zoneType);
    tecZoneGetIJK(infH, 1, &iMax, &jMax, &kMax);
    if ((zoneType != 5) && (zoneType != 3))
    {
        printf("Zone is weird/Not supported\n");
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
    tecFileWriterOpen((fpref + to_string(fileID) + suffix).c_str(), "Code out", varstr.c_str(), 1, 0, 1, NULL, &outfH);
    assert(outfH != NULL);
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

    if (GEMSbin)
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

void tecIO::addVar(string var, string &norm)
{
    varName.push_back(var);
    int varI = getVariableIndex(var, prefix + to_string(snap0) + suffix);
    varIndex.push_back(varI);
    normID.push_back(norm);
    double temp;
    temp = stod(norm);
    scalingInput.push_back(temp);
    assert(varName.size() == varIndex.size());
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
    assert(fH != NULL);
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

void tecIO::genHash(string filename)
{
    if (idx.size() != 0)
        return;

    idx.resize(nCells, 0);
    iota(idx.begin(), idx.end(), 0);
    vector<int> cellID(nCells);
    void *fH;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (filename != "")
    {
        int var_index = 0;
        var_index = getVariableIndex("cell_id", filename);

        // faster if every process loads cell_id instead of broadcasting
        int hashType;
        tecFileReaderOpen(filename.c_str(), &fH);
        tecZoneVarGetType(fH, 1, var_index, &hashType);
        if (hashType == 3)
        {
            tecZoneVarGetInt32Values(fH, 1, var_index, 1, nCells, cellID.data());
        }
        else
        {
            MPI_Abort(MPI_COMM_WORLD,-1);
        }
        tecFileReaderClose(&fH);
        for (int i = 0; i < nCells; i++)
        {
            cellID[i]--;
        }

        cout << "Hash table built" << endl;
        stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return cellID[i] < cellID[j]; });
        cout << "Finish sort" << endl;
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
    genHash(prefix + to_string(snap0) + suffix);
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
            cout << "Invalid centering file: " << centerMethod << endl;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        cout << "Centering calculation has not been re-implemented" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // update global param
    centeringIsField = isField;

    // expand constants to full field for convenience sake
    if (!isField)
    {
        cout << "Distributing scalar centering values has not been re-implemented" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // write centering field to file
    cout << "Centering calculated" << endl;
    if (writeToDisk)
    {
        // TODO: get SZPLT to output correctly, without silly fileID requirement
        centerVec->write_ascii("centerProf.dat", "centerProf");
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
        cout << "Centering isn't setup, call calcCentering first" << endl;
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
    genHash(prefix + to_string(snap0) + suffix);
    scalingSubVec = new pMat(nPoints, 1, dataMat->pG, 0, 0, 0.0, false);
    scalingDivVec = new pMat(nPoints, 1, dataMat->pG, 0, 0, 0.0, false);
    scalingSubVecFull = new pMat(nPoints, 1, dataMat->pG, 0, 0, 0.0, false);

    if ((scaleMethod == "sphere") && (isField))
    {
        cout << "If choosing scaleMethod = sphere, you must set scaleIsField = false" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
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
            cout << "Could not find scaling files with prefix  " << scaleMethod << endl;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        if (calc)
        {
            cout << "Scaling calculation has not been re-implemented" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
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
        scalingSubVec->write_ascii("scalingSubProf.dat", "scalingSubProf");
        scalingDivVec->write_ascii("scalingDivProf.dat", "scalingDivProf");

        // if also centering, combine into single profile (useful for GEMS)
        if (isCentered)
        {
            scalingSubVecFull->matrix_Sum('N', nPoints, 1, centerVec, 0, 0, 0, 0, 1.0, 1.0);
            scalingSubVecFull->matrix_Sum('N', nPoints, 1, scalingSubVec, 0, 0, 0, 0, 1.0, 1.0);
            cout << "Writing combined subtractive scaling factors" << endl;
            scalingSubVecFull->write_ascii("scalingSubProfFull.dat", "scalingSubProfFull");
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
        cout << "Scaling isn't setup, call calcScaling first" << endl;
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

void tecIO::activateGEMSbin(string file)
{
    GEMSbin = true;
    genHash(file);
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
