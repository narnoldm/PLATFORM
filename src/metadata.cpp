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
    fread(&(header[0]), sizeof(int), 1, fid);
    fread(&(header[1]), sizeof(int), 1, fid);
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

bool meta::readSingle(int fileID, double *point)
{
    cout << "meta read " << fileID << endl;
    cout<< (prefix + to_string(fileID) + suffix)<<endl;
    FILE *fid;
    fid = fopen((prefix + to_string(fileID) + suffix).c_str(), "rb");
    int header[2] = {0, 0};
    fread(&(header[0]), sizeof(int), 1, fid);
    fread(&(header[1]), sizeof(int), 1, fid);
    assert(header[1] == 1);
    assert(header[0] == nPoints);
    fread(point, sizeof(double), nPoints, fid);
    fclose(fid);
    return true;
}

bool meta::batchRead(pMat *loadMat)
{
    double t1,t2; 
    t1=MPI_Wtime();
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

                cout << "proc " << iP << " is reading file " << fileIndex << endl;

                readSingle(fileIndex, loadMat->dataD.data() + nPoints * localC);
                localC++;
            }
        }
        miscProcessing(loadMat);
        cout << "waiting for other processes read" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        cout << "even Read" << endl;
        int currentCol = 0;
        vector<double> tempR;
        tempR.resize(nPoints);
        for (int j = 0; j < nSets; j++)
        {
            if (loadMat->pG->mycol == (j / loadMat->nb) % loadMat->pG->pcol)
            {
                readSingle(snap0 + j * snapSkip, tempR.data());
                for (int i = 0; i < nPoints; i++)
                {
                    int xi = i % loadMat->mb;
                    int li = i / (loadMat->pG->prow * loadMat->mb);
                    if (loadMat->pG->myrow == (i / loadMat->mb) % loadMat->pG->prow)
                    {
                        loadMat->dataD[currentCol * loadMat->myRC[0] + xi + li * loadMat->mb] = tempR[i];
                    }
                }
                currentCol++;
            }
        }
        tempR.clear();
    }
    t2=MPI_Wtime();
    cout<<"batch Read took "<<t2-t1<<" secs"<<endl;
}

bool meta::batchRead(pMat *loadMat, int ii)
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

            cout << "proc " << iP << " is reading file " << fileIndex << "\r";

            readSingle(fileIndex, loadMat->dataD.data() + nPoints * localC);
            localC++;
        }
        miscProcessing(loadMat);
        cout << "waiting for other processes read" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else
    {
        cout << "even Read single" << endl;
        int currentCol = 0;
        vector<double> tempR;
        tempR.resize(nPoints);
        int j = ii;
        if (loadMat->pG->mycol == (j / loadMat->nb) % loadMat->pG->pcol)
        {
            readSingle(snap0 + j * snapSkip, tempR.data());
            for (int i = 0; i < nPoints; i++)
            {
                int xi = i % loadMat->mb;
                int li = i / (loadMat->pG->prow * loadMat->mb);
                if (loadMat->pG->myrow == (i / loadMat->mb) % loadMat->pG->prow)
                {
                    loadMat->dataD[currentCol * loadMat->myRC[0] + xi + li * loadMat->mb] = tempR[i];
                }
            }
            currentCol++;
        }
        tempR.clear();
    }
}

bool meta::writeSingle(int fileID, double *point, string fpref)
{
    cout << "meta write single" << endl;
    FILE *fid;
    fid = fopen((fpref + to_string(fileID) + suffix).c_str(), "wb");

    const int ONE = 1;
    int points = nPoints;
    fwrite(&points, sizeof(int), 1, fid);
    fwrite(&ONE, sizeof(int), 1, fid);
    fwrite(point, sizeof(double), nPoints, fid);
    fclose(fid);
    return true;
}

bool meta::batchWrite(pMat *loadMat)
{
    batchWrite(loadMat, "out/", prefix);
}
bool meta::batchWrite(pMat *loadMat, string dir, string fpref)
{
    batchWrite(loadMat, dir, fpref, 0, nSets, 1);
}
bool meta::batchWrite(pMat *loadMat, string dir, string fpref, int nModes)
{
    batchWrite(loadMat, dir, fpref, 0, nModes, 1);
}
bool meta::batchWrite(pMat *loadMat, string dir, string fpref, int mStart, int mEnd, int mSkip)
{
    batchWrite(loadMat, dir, fpref, mStart, mEnd, mSkip, snap0, snapSkip);
}

bool meta::batchWrite(pMat *loadMat, string dir, string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip) {
    
    assert(system(NULL)); //check if system commands work
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
        system(("mkdir " + dir).c_str());

    MPI_Barrier(MPI_COMM_WORLD);
    int iP = 0, fileIndex, localC = 0;

    if (!isInit) {
        nPoints = loadMat->M;
        nSets = loadMat->N;
    }

    double t1, t2;
    t1 = MPI_Wtime();
    if (loadMat->block == 1)
    {
        assert(loadMat->mb == nPoints);
        for (int i = mStart; i < mEnd; i = i + mSkip)
        {
            iP = (int)(i / loadMat->nb);
            while (iP > (loadMat->pG->pcol - 1))
            {
                iP = iP - loadMat->pG->pcol;
            }
            if (loadMat->pG->rank == iP)
            {
                fileIndex = fStart + (i - mStart) * fSkip;

                cout << "proc " << iP << " is writing file " << fileIndex << "\r";
                writeSingle(fileIndex, loadMat->dataD.data() + nPoints * localC, dir + "/" + fpref);
                localC++;
            }
        }
    }
    else
    {
        int currentCol = 0;
        vector<double> tempR;
        tempR.resize(nPoints, 0);
        MPI_Comm col_comms;
        loadMat->commCreate(col_comms, 0);
        for (int j = mStart; j < mEnd; j = j + mSkip)
        {
            std::fill(tempR.begin(), tempR.end(), 0.0);
            if (loadMat->pG->mycol == (j / loadMat->nb) % loadMat->pG->pcol)
            {

                for (int i = 0; i < nPoints; i++)
                {
                    int xi = i % loadMat->mb;
                    int li = i / (loadMat->pG->prow * loadMat->mb);
                    if (loadMat->pG->myrow == (i / loadMat->mb) % loadMat->pG->prow)
                    {
                        tempR[i] = loadMat->dataD[currentCol * loadMat->myRC[0] + xi + li * loadMat->mb];
                    }
                }
                currentCol++;
            }
            MPI_Allreduce(MPI_IN_PLACE, tempR.data(), tempR.size(), MPI_DOUBLE, MPI_SUM, col_comms);
            if ((loadMat->pG->mycol == (j / loadMat->nb) % loadMat->pG->pcol) && (loadMat->pG->myrow == 0))
            {
                fileIndex = fStart + (j - mStart) * fSkip;
                printf("proc %d is writing %d\n", loadMat->pG->rank, fileIndex);
                writeSingle(fileIndex, tempR.data(), dir + "/" + fpref);
            }
        }
        tempR.clear();
        cout << endl;
        MPI_Comm_free(&col_comms);
    }
    t2 = MPI_Wtime();
    cout << "batch Write took " << t2 - t1 << " secs" << endl;
}

void meta::miscProcessing(pMat *Mat)
{
    cout << "no additional processing for binary" << endl;
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
            tecFileReaderOpen((prefix + std::to_string(i) + suffix).c_str(), &fH);
            assert(fH != NULL);
            tecFileReaderClose(&fH);
        }
    }
}
bool tecIO::readSingle(int fileID, double *point)
{
    cout << (prefix + std::to_string(fileID) + suffix) << "\r";
    void *fH = NULL;
    tecFileReaderOpen((prefix + std::to_string(fileID) + suffix).c_str(), &fH);
    int type;
    assert(fH != NULL);
    std::vector<float> get;
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
            get.clear();
        }
        else if (type == 2)
        {
            tecZoneVarGetDoubleValues(fH, 1, varIndex[i], 1, nCells, &(point[i * nCells]));
        }
        if (reorder)
        {
            //cout << "reording slice" << endl;
            std::vector<double> temp(nCells, 0.0);
            for (int j = 0; j < nCells; j++)
            {
                temp[j] = point[i * nCells + j];
            }
            for (int j = 0; j < nCells; j++)
            {
                point[i * nCells + j] = temp[idx[j]];
                //hash[j]=j;
            }
            temp.clear();
        }
    }
    tecFileReaderClose(&fH);
}

bool tecIO::writeSingle(int fileID, double *point, string fpref)
{
    void *infH = NULL;
    void *outfH = NULL;
    if (fixedMesh)
    {
        tecFileReaderOpen((meshFile).c_str(), &infH);
    }
    else
    {
        tecFileReaderOpen((prefix + std::to_string(fileID) + suffix).c_str(), &infH);
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

    std::string varstr = "";
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
    tecFileWriterOpen((fpref + std::to_string(fileID) + suffix).c_str(), "Code out", varstr.c_str(), 1, 0, 1, NULL, &outfH);
    assert(outfH != NULL);
    std::vector<int> varTypes(dim + numVars);
    std::vector<int> valueLoc(dim + numVars);
    std::vector<int> passive(dim + numVars);
    std::vector<int> shareVar(dim + numVars);
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
    tecZoneCreateFE(outfH, std::to_string(fileID).c_str(), zoneType, iMax, jMax, &varTypes[0], &shareVar[0], &valueLoc[0], &passive[0], shareCon, 0, 0, &outZone);
    tecZoneSetUnsteadyOptions(outfH, outZone, fileID, (int)fileID);

    vector<float> nfDat;
    vector<double> ndDat;
    for (int i = 0; i < dim; i++)
    {
        if (varTypes[i] == 1) //float
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
    nfDat.clear();
    ndDat.clear();
    for (int i = dim; i < (dim + numVars); i++)
    {
        if (reorder)
        {
            vector<double> temp(jMax, 0.0);
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
    long numValues;
    tecZoneNodeMapGetNumValues(infH, 1, jMax, &numValues);
    vector<int> nodeMap(numValues);
    tecZoneNodeMapGet(infH, 1, 1, jMax, nodeMap.data());
    tecZoneNodeMapWrite32(outfH, 1, 0, 1, numValues, nodeMap.data());
    nodeMap.clear();
    tecFileReaderClose(&infH);
    tecFileWriterClose(&outfH);
    varTypes.clear();
    valueLoc.clear();
    passive.clear();
    shareVar.clear();

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
void tecIO::miscProcessing(pMat *Mat)
{
}
bool tecIO::writeSingleFile(std::string filename, std::vector<std::string> &fvars, double *point, std::string meshPfile)
{
    void *infH = NULL;
    void *outfH = NULL;
    tecFileReaderOpen((meshPfile).c_str(), &infH);

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

    std::string varstr = "";
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
    for (int i = 0; i < fvars.size(); i++)
    {
        varstr = varstr + "," + fvars[i];
    }
    tecFileWriterOpen(filename.c_str(), "Code out", varstr.c_str(), 1, 0, 1, NULL, &outfH);
    assert(outfH != NULL);
    std::vector<int> varTypes(dim + fvars.size());
    std::vector<int> valueLoc(dim + fvars.size());
    std::vector<int> passive(dim + fvars.size());
    std::vector<int> shareVar(dim + fvars.size());
    for (int i = 0; i < (dim + fvars.size()); i++)
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
    tecZoneCreateFE(outfH, std::to_string(snap0).c_str(), zoneType, iMax, jMax, &varTypes[0], &shareVar[0], &valueLoc[0], &passive[0], shareCon, 0, 0, &outZone);
    tecZoneSetUnsteadyOptions(outfH, outZone, snap0, (int)snap0);

    vector<float> nfDat;
    vector<double> ndDat;
    for (int i = 0; i < dim; i++)
    {
        if (varTypes[i] == 1) //float
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
    nfDat.clear();
    ndDat.clear();
    for (int i = dim; i < (dim + fvars.size()); i++)
    {
        if (reorder)
        {
            vector<double> temp(jMax, 0.0);
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
    long numValues;
    tecZoneNodeMapGetNumValues(infH, 1, jMax, &numValues);
    vector<int> nodeMap(numValues);
    tecZoneNodeMapGet(infH, 1, 1, jMax, nodeMap.data());
    tecZoneNodeMapWrite32(outfH, 1, 0, 1, numValues, nodeMap.data());
    nodeMap.clear();
    tecFileReaderClose(&infH);
    tecFileWriterClose(&outfH);
    varTypes.clear();
    valueLoc.clear();
    passive.clear();
    shareVar.clear();
}
void tecIO::addVar(string var, string &norm)
{
    varName.push_back(var);
    int varI = getVariableIndex(var, prefix + to_string(snap0) + suffix);
    varIndex.push_back(varI);
    normID.push_back(norm);
    double temp;
    temp = stod(norm);
    normFactor.push_back(temp);
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
            checkMeshDim((prefix + std::to_string(snap0) + suffix));
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
    cellID.resize(nCells, 0);
    void *fH;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (filename != "")
    {
        int var_index = 0;
        var_index = getVariableIndex("cell_id", filename);
        if (!rank)
        {
            int hashType;
            cellID.resize(nCells);
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

            printf("hash table built\nDistributing to procs\n");
        }

        MPI_Bcast(idx.data(), nCells, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(cellID.data(), nCells, MPI_INT, 0, MPI_COMM_WORLD);
        stable_sort(idx.begin(), idx.end(), [&](int i, int j) { return cellID[i] < cellID[j]; });
    }
    else
    {
        cout << "assume default hash" << endl;
        for (int i = 0; i < nCells; i++)
        {
            cellID[i] = i;
        }
    }
}

void tecIO::normalize(pMat *dataMat)
{
    cout << "Normalizing Matrix by Norm Factor" << endl;
    if (dataMat->block == 1)
    {
        int numFiles = dataMat->nelements / nPoints;
        cout << "proc " << dataMat->pG->rank << "has " << numFiles << " Files " << endl;
        for (int i = 0; i < dataMat->N; i++)
        {
            for (int j = 0; j < numVars; j++)
            {
                for (int k = 0; k < nCells; k++)
                {
                    dataMat->dataD[i * nPoints + j * nCells + k] /= normFactor[j];
                }
            }
        }
    }
    else
    {
        int currentCol = 0;
        for (int j = 0; j < dataMat->N; j++)
        {
            if (dataMat->pG->mycol == (j / dataMat->nb) % dataMat->pG->pcol)
            {
                for (int i = 0; i < nPoints; i++)
                {
                    int xi = i % dataMat->mb;
                    int li = i / (dataMat->pG->prow * dataMat->mb);
                    if (dataMat->pG->myrow == (i / dataMat->mb) % dataMat->pG->prow)
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] /= normFactor[i / nCells];
                    }
                }
                currentCol++;
            }
        }
    }
}

void tecIO::unNormalize(pMat *dataMat)
{
    cout << "UNNormalizing Matrix by Norm Factor" << endl;
    if (dataMat->block == 1)
    {
        int numFiles = dataMat->nelements / nPoints;
        cout << "proc " << dataMat->pG->rank << "has " << numFiles << " Files " << endl;
        for (int i = 0; i < dataMat->N; i++)
        {
            for (int j = 0; j < numVars; j++)
            {
                for (int k = 0; k < nCells; k++)
                {
                    dataMat->dataD[i * nPoints + j * nCells + k] *= normFactor[j];
                }
            }
        }
    }
    else
    {
        int currentCol = 0;
        for (int j = 0; j < dataMat->N; j++)
        {
            if (dataMat->pG->mycol == (j / dataMat->nb) % dataMat->pG->pcol)
            {
                for (int i = 0; i < nPoints; i++)
                {
                    int xi = i % dataMat->mb;
                    int li = i / (dataMat->pG->prow * dataMat->mb);
                    if (dataMat->pG->myrow == (i / dataMat->mb) % dataMat->pG->prow)
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] *= normFactor[i / nCells];
                    }
                }
                currentCol++;
            }
        }
    }
}

void tecIO::calcNorm(pMat *dataMat)
{
    cout << "calculating norms" << endl;
    if (dataMat->block == 1)
    {
        int numFiles = dataMat->nelements / nPoints;
        std::vector<double> mag(nCells * numFiles, 0.0);
        double maxmag = 0.0;
        double group = 0.0;
        for (int i = 0; i < numVars; i++)
        {
            if (normFactor[i] < 0)
            {
                group = normFactor[i];
                for (int m = 0; m < (nCells * numFiles); m++)
                    mag[m] = 0.0;
                maxmag = 0.0;

                cout << varName[i] << " is in normalization group " << normFactor[i] << endl;
                for (int j = 0; j < numVars; j++)
                {
                    if (normFactor[j] == group)
                    {
                        for (int f = 0; f < numFiles; f++)
                        {
                            for (int n = 0; n < nCells; n++)
                            {
                                mag[f * nCells + n] += (dataMat->dataD[f * numVars * nCells + j * nCells + n]) * (dataMat->dataD[f * numVars * nCells + j * nCells + n]);
                            }
                        }
                    }
                }

                for (int n = 0; n < (nCells * numFiles); n++)
                {
                    mag[n] = sqrt(mag[n]);
                    if (mag[n] > maxmag)
                        maxmag = mag[n];
                }
                for (int j = 0; j < numVars; j++)
                {
                    if (normFactor[j] == group)
                    {
                        normFactor[j] = maxmag;
                        cout << "Setting norm factor " << j << " to " << maxmag << endl;
                    }
                }
            }
            else
            {
                cout << varName[i] << " is using specified normalization of " << normFactor[i] << endl;
            }
        }
        mag.clear();
        for (int i = 0; i < numVars; i++)
        {
            MPI_Allreduce(MPI_IN_PLACE, normFactor.data() + i, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            cout << "Synched norm factor for " << varName[i] << " is : " << setprecision(numeric_limits<double>::digits10) << normFactor[i] << endl;
        }
        if(dataMat->pG->rank == 0)
        {
                std::vector<double> normVec(numVars * nCells);
                for (int i = 0; i < numVars; i++)
                {
                        for (int j = 0; j < nCells; j++)
                        {
                                normVec[i * nCells + j] = normFactor[i];
                        }
                }
                writeSingle(snap0, normVec.data(), "norm");
                printASCIIVecP0("norm.data",normVec,normVec.size());
                normVec.clear();
                //delete [] normVec;
        }
    }
    else
    {

        for (int k = 0; k < numVars; k++)
        {
            if (normFactor[k] < 0)
            {
                double maxmag = 0;
                for (int j = 0; j < nSets; j++)
                {

                    vector<double> mag(nCells, 0.0);
                    std::fill(mag.begin(), mag.end(), 0.0);
                    for (int i = 0; i < nCells; i++)
                    {
                        for (int g = 0; g < numVars; g++)
                        {

                            if (normFactor[g] == normFactor[k])
                            {

                                mag[i] += std::pow(dataMat->getLocalElement(g * nCells + i, j), 2);
                            }
                        }
                    }
                    MPI_Allreduce(MPI_IN_PLACE, mag.data(), nCells, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
                    for (int i = 0; i < nCells; i++)
                    {
                        if (mag[i] > maxmag)
                        {
                            maxmag = mag[i];
                        }
                    }
                }

                maxmag = sqrt(maxmag);
                for (int g = 0; g < numVars; g++)
                {
                    if (g == k)
                    {
                        continue;
                    }
                    else if (normFactor[g] == normFactor[k])
                    {
                        normFactor[g] = maxmag;
                    }
                }
                normFactor[k] = maxmag;
            }
            cout << "Synched norm factor for " << varName[k] << " is : " << setprecision(numeric_limits<double>::digits10) << normFactor[k] << endl;
        }
        /*if(dataMat->pG->rank == 0)
        {
                std::vector<double> normVec(numVars * nCells);
                for (int i = 0; i < numVars; i++)
                {
                        for (int j = 0; j < nCells; j++)
                        {
                                normVec[i * nCells + j] = normFactor[i];
                        }
                }
                writeSingle(snap0, normVec.data(), "norm");
                printASCIIVecP0("norm.data",normVec.data(),normVec.size());
                normVec.clear();
        }*/
    }
}

void tecIO::subAvg(pMat *dataMat)
{
    if (average.size() == 0)
    {
        cout << "Average isn't setup call calcAvg first" << endl;
        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    cout << "Subtracting Average" << endl;
    if (dataMat->block == 1)
    {
        int numFiles = dataMat->nelements / nPoints;
        for (int i = 0; i < dataMat->N; i++)
        {
            for (int j = 0; j < nPoints; j++)
                dataMat->dataD[i * nPoints + j] -= average[j];
        }
    }
    else
    {
        int currentCol = 0;
        for (int j = 0; j < dataMat->N; j++)
        {
            if (dataMat->pG->mycol == (j / dataMat->nb) % dataMat->pG->pcol)
            {
                for (int i = 0; i < nPoints; i++)
                {
                    int xi = i % dataMat->mb;
                    int li = i / (dataMat->pG->prow * dataMat->mb);
                    if (dataMat->pG->myrow == (i / dataMat->mb) % dataMat->pG->prow)
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] -= average[i];
                    }
                }
                currentCol++;
            }
        }
    }
    cout << "average subtracted" << endl;
}

void tecIO::addAvg(pMat *dataMat)
{
    if (average.size() == 0)
    {
        cout << "Average isn't setup call calcAvg first" << endl;
        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    cout << "Adding Average" << endl;
    if (dataMat->block == 1)
    {
        int numFiles = dataMat->nelements / nPoints;
        for (int i = 0; i < dataMat->N; i++)
        {
            for (int j = 0; j < nPoints; j++)
                dataMat->dataD[i * nPoints + j] += average[j];
        }
    }
    else
    {
        int currentCol = 0;
        for (int j = 0; j < dataMat->N; j++)
        {
            if (dataMat->pG->mycol == (j / dataMat->nb) % dataMat->pG->pcol)
            {
                for (int i = 0; i < nPoints; i++)
                {
                    int xi = i % dataMat->mb;
                    int li = i / (dataMat->pG->prow * dataMat->mb);
                    if (dataMat->pG->myrow == (i / dataMat->mb) % dataMat->pG->prow)
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] += average[i];
                    }
                }
                currentCol++;
            }
        }
    }
    cout << "average added" << endl;
}

void tecIO::calcAvg(pMat *dataMat)
{

    if (average.size() != nPoints)
    {
        cout << "Allocating Average as " << nPoints << " data Points" << endl;
        average.resize(nPoints, 0.0);
        cout << "Average Allocated" << endl;
    }
    std::fill(average.begin(), average.end(), 0.0);

    if (dataMat->block == 1)
    {
        int numFiles = dataMat->nelements / nPoints;
        for (int i = 0; i < nPoints; i++)
        {
            for (int j = 0; j < numFiles; j++)
            {
                average[i] += dataMat->dataD[j * nPoints + i];
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, average.data(), nPoints, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int i = 0; i < nPoints; i++)
            average[i] /= nSets;

        if (dataMat->pG->rank == 0)
            writeSingle(snap0, average.data(), "average1");
    }
    else
    {
        for (int i = 0; i < nPoints; i++)
        {
            for (int j = 0; j < nSets; j++)
            {

                average[i] += dataMat->getLocalElement(i, j);
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, average.data(), nPoints, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        for (int i = 0; i < nPoints; i++)
        {
            average[i] /= nSets;
        }
        cout << "average calculated" << endl;
        if (dataMat->pG->rank == 0)
            writeSingle(snap0, average.data(), "average0");
        cout << "average outputed" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void tecIO::readAvg(std::string filename)
{
    if (average.size() != nPoints)
    {
        cout << "Allocating Average as " << nPoints << " cells" << endl;
        average.resize(nPoints, 0.0);
        cout << "Average Allocated" << endl;
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //if (rank == 0)
    //{
    void *fH;
    tecFileReaderOpen(filename.c_str(), &fH);
    int type;
    std::vector<float> get;
    int ii;
    for (int i = 0; i < numVars; i++)
    {
        ii = getVariableIndex(varName[i], filename);
        tecZoneVarGetType(fH, 1, ii, &type);
        if (type == 1)
        {
            //get= new float[nCells];
            get.resize(nCells);
            tecZoneVarGetFloatValues(fH, 1, ii, 1, nCells, get.data());
            for (int j = 0; j < nCells; j++)
            {
                average[j + i * nCells] = (double)get[j];
            }
            //delete [] get;
            get.clear();
        }
        else if (type == 2)
        {
            tecZoneVarGetDoubleValues(fH, 1, ii, 1, nCells, &(average[i * nCells]));
        }
        if (reorder)
        {
            //cout << "reording slice" << endl;
            std::vector<double> temp(nCells, 0.0);
            for (int j = 0; j < nCells; j++)
            {
                temp[j] = average[i * nCells + j];
            }
            for (int j = 0; j < nCells; j++)
            {
                average[i * nCells + j] = temp[idx[j]];
                //hash[j]=j;
            }
            temp.clear();
        }
    }
    tecFileReaderClose(&fH);
    std::cout << "average loaded from :" << filename << std::endl;
    //}
    //printf("%d\n",nPoints);
    genHash(filename);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cout << "outputing ascii centering" << std::endl;
    if (!rank)
    {
        FILE *fid;
        string asciiName = "centerProf.dat";
        if ((fid = fopen(asciiName.c_str(), "w")) == NULL)
        {
            printf("error with file open\n");
        }
        fprintf(fid, "file= %s\n", asciiName.c_str());
        if (reorder)
        {
            for (int i = 0; i < average.size(); i++)
            {
                fprintf(fid, "%16.16E\n", average[i]);
            }
        }
        else
        {
            for (int i = 0; i < numVars; i++)
            {
                for (int j = 0; j < nCells; j++)
                {
                    fprintf(fid, "%16.16E\n", average[i * nCells + idx[j]]);
                }
            }
        }
    }
    std::cout << "average Broad-casted" << std::endl;
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
int compareMeta(meta* meta1, meta* meta2) {

	// check path members
	if ( (meta1->prefix == meta2->prefix) && (meta1->suffix == meta2->suffix) ) {

		if ( (meta1->snap0 == meta2->snap0) && 
			 (meta1->snapF == meta2->snapF) &&
			 (meta1->snapSkip == meta2->snapSkip) ) {
			
			// objects are identical
			return(1);

		} else {
			// objects differ only in the file counts
			return(2);

		}
		return(true);

	} else {
		// objects are totally different
		return(0);
	}

}
