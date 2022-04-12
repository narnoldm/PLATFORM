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

bool meta::batchWrite(pMat *loadMat, string dir, string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip)
{

    assert(system(NULL)); //check if system commands work
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
        system(("mkdir " + dir).c_str());

    MPI_Barrier(MPI_COMM_WORLD);
    int iP = 0, fileIndex, localC = 0;

    if (!isInit)
    {
        nPoints = loadMat->M;
        nSets = loadMat->N;
    }

    double t1, t2;
    t1 = MPI_Wtime();
    int currentCol = 0;
    vector<double> tempR;
    tempR.resize(nPoints, 0);
    MPI_Comm col_comms;
    loadMat->commCreate(col_comms, 0);
    for (int j = mStart; j < mEnd; j = j + mSkip)
    {
        fill(tempR.begin(), tempR.end(), 0.0);
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

    isCentered = true;
    genHash(prefix + to_string(snap0) + suffix);

    // update this check when a new method is added
    if ((centerMethod != "avg") && (centerMethod != "avgmag"))
    {

        // load from file
        // TODO: allow reading from small vector
        // TODO: allow reading from ASCII file
        isField = true;
        if (centerMethod.substr(centerMethod.size()-6, 6) == ".szplt")
        {
            readSZPLTToVec(centerMethod, centerVec);
        }
        else if (centerMethod.substr(centerMethod.size()-4, 4) == ".dat")
        {
            readDATToVec(centerMethod, centerVec);
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
        // otherwise, calculate centering profile
        int centerSize;
        if (isField)
        {
            centerSize = nPoints;
        }
        else
        {
            centerSize = numVars;
        }
        cout << "Allocating centering as " << centerSize << " values" << endl;
        centerVec.resize(centerSize, 0.0);
        cout << "Centering allocated" << endl;

        double val;
        vector<double> valVec(nCells);
        bool skip;
        vector<int> skipFlags;
        for (int k = 0; k < numVars; ++k)
        {
            if (scalingInput[k] < 0)
            {
                // check if this group has already been evaluated
                skip = false;
                for (int g = 0; g < skipFlags.size(); ++g)
                {
                    if (scalingInput[k] == skipFlags[g])
                    {
                        skip = true;
                        break;
                    }
                }
                if (skip)
                {
                    continue;
                }

                if (centerMethod == "avg")
                {
                    calcGroupQuant(dataMat, val, valVec, k, "avg", isField);
                }
                else if (centerMethod == "avgmag")
                {
                    calcGroupQuant(dataMat, val, valVec, k, "avgmag", isField);
                }

                // final calculations for centering method
                if ((centerMethod == "avg") || (centerMethod == "avgmag"))
                {
                    if (isField)
                    {
                        for (int i = 0; i < nCells; ++i)
                        {
                            centerVec[k * nCells + i] = valVec[i];
                        }
                    }
                    else
                    {
                        centerVec[k] = val;
                    }
                }

                // distribute to fields in same group instead of recalculating
                for (int g = 0; g < numVars; g++)
                {
                    if (g == k)
                    {
                        continue;
                    }
                    else if (scalingInput[g] == scalingInput[k])
                    {
                        if (isField)
                        {
                            for (int i = 0; i < nCells; ++i)
                            {
                                centerVec[g * nCells + i] = centerVec[k * nCells + i];
                            }
                        }
                        else
                        {
                            centerVec[g] = centerVec[k];
                        }
                        skipFlags.push_back(k);
                    }
                }
            }

            if (!isField)
            {
                cout << "Centering factor for " << varName[k] << " is: " << setprecision(numeric_limits<double>::digits10) << centerVec[k] << endl;
            }
            else
            {
                cout << "Centering field calculation for " << varName[k] << " complete" << endl;
            }
        }
    }

    // update global param
    centeringIsField = isField;

    // expand constants to full field for convenience sake
    if (!isField)
    {
        vector<double> centerVals(numVars);
        copy(centerVec.begin(), centerVec.end(), centerVals.begin());
        centerVec.resize(nPoints, 0.0);
        for (int k = 0; k < numVars; ++k) 
        {
            for (int i = 0; i < nCells; ++i)
            {
                centerVec[k * nCells + i] = centerVals[k];
            }
        }
    }

    // write centering field to file
    cout << "Centering calculated" << endl;
    if (dataMat->pG->rank == 0)
    {
        // TODO: get SZPLT to output correctly, without silly fileID requirement
        // writeSingle(0, centerVec.data(), "centerProf");
        // convert to cell_id order and write
        vector<double> vecOut(nPoints, 0.0);
        vecToCellIDOrder(centerVec, vecOut);
        writeASCIIDoubleVec("centerProf.dat", vecOut);
    }
    cout << "Centering files written" << endl;
    MPI_Barrier(MPI_COMM_WORLD);

}

void tecIO::centerData(pMat *dataMat)
{
    centerData(dataMat, false);
}

void tecIO::centerData(pMat *dataMat, bool uncenter)
{
    if (centerVec.size() == 0)
    {
        cout << "Centering isn't setup, call calcCentering first" << endl;
        MPI_Abort(MPI_COMM_WORLD,-1);
    }

    if (uncenter)
    {
        cout << "Uncentering" << endl;
    }
    else
    {
        cout << "Centering" << endl;
    }

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
                    // center data
                    if (uncenter)
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] += centerVec[i];
                    }
                    else
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] -= centerVec[i];
                    }
                }
            }
            currentCol++;
        }
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

void tecIO::calcScaling(pMat *dataMat, string scaleMethod, bool isField)
{

    isScaled = true;

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
        if (scalingInput[i] > 0)
        {
            if (!isField)
            {
                cout << "WARNING: Passed a scaling factor, all scaling will by by scalars" << endl;
            }
            isField = false;
            break;
        }
    }

    // update this check when a new method is added
    if ((scaleMethod != "minmax") &&
        (scaleMethod != "standardize") &&
        (scaleMethod != "sphere"))
    {
        // load from file
        // TODO: allow reading from small vector
        // TODO: allow reading from ASCII file
        isField = true;
        ifstream szlFile((scaleMethod + "SubProf.szplt").c_str());
        ifstream datFile((scaleMethod + "SubProf.dat").c_str());
        if (szlFile.good())
        {
            readSZPLTToVec(scaleMethod + "SubProf.szplt", scalingSubVec);
            readSZPLTToVec(scaleMethod + "DivProf.szplt", scalingDivVec);
        }
        else if (datFile.good())
        {
            readDATToVec(scaleMethod + "SubProf.dat", scalingSubVec);
            readDATToVec(scaleMethod + "DivProf.dat", scalingDivVec);
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

        // initialize scaling vectors to appropriate size
        int scaleVecSize;
        if (isField)
        {
            scaleVecSize = nPoints;
        }
        else
        {
            scaleVecSize = numVars;
        }
        cout << "Allocating scaling as " << scaleVecSize << " values" << endl;
        scalingSubVec.resize(scaleVecSize, 0.0);
        scalingDivVec.resize(scaleVecSize, 1.0);
        cout << "Scaling allocated" << endl;

        // need a copy of data for standardization
        pMat* dataMatCopy;
        if (scaleMethod == "standardize")
        {
            dataMatCopy = new pMat(nPoints, nSets, dataMat->pG, 0, 0, 0.0, false);
        }

        double val1, val2;
        vector<double> valVec1(nCells);
        vector<double> valVec2(nCells);
        bool skip;
        vector<int> skipFlags;
        for (int k = 0; k < numVars; ++k)
        {
            if (scalingInput[k] < 0)
            {
                // check if this group has already been evaluated
                skip = false;
                for (int g = 0; g < skipFlags.size(); ++g)
                {
                    if (scalingInput[k] == skipFlags[g])
                    {
                        skip = true;
                        break;
                    }
                }
                if (skip)
                {
                    continue;
                }

                if (scaleMethod == "minmax")
                {
                    calcGroupQuant(dataMat, val1, valVec1, k, "min", isField);
                    calcGroupQuant(dataMat, val2, valVec2, k, "max", isField);
                }
                else if (scaleMethod == "standardize")
                {
                    // requires a few extra steps

                    // calculate average
                    calcGroupQuant(dataMat, val1, valVec1, k, "avg", isField);
                    pMat* avgVecP0 = new pMat(nCells, 1, dataMat->pG, 0, 2, 0.0, false);
                    pMat* avgVec = new pMat(nCells, 1, dataMat->pG, 0, 0, 0.0, false);
                    int rank;
                    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                    if (!rank)
                    {
                        for (int i = 0; i < nCells; ++i)
                        {
                            if (isField)
                            {
                                avgVecP0->dataD[i] = valVec1[i];
                            }
                            else
                            {
                                avgVecP0->dataD[i] = val1;
                            }
                        }
                    }
                    avgVec->changeContext(avgVecP0, false);

                    // subtract average
                    // TODO: the repeated changeContext and squaring of all elements is CRAZY inefficient
                    dataMatCopy->changeContext(dataMat, false);
                    for (int j = 0; j < nSets; ++j)
                    {
                        for (int g = 0; g < numVars; ++g)
                        {
                            if (scalingInput[k] == scalingInput[g])
                            {
                                dataMatCopy->matrix_Sum('N', nCells, 1, avgVec, 0, 0, g * nCells, j, -1.0, 1.0);
                            }
                        }
                    }

                    // TODO: this is subtly incorrect for grouped variables
                    // This computes the variance as sum((x1 - mu)^2 + (x2 - mu)^2) / N
                    // Should be sum((x1 + x2 - 2 * mu)^2) / N
                    // Not sure how this could be done efficiently

                    // square elements
                    for (int i = 0; i < dataMatCopy->dataD.size(); ++i)
                    {
                        dataMatCopy->dataD[i] = pow(dataMatCopy->dataD[i], 2);
                    }
                    MPI_Barrier(MPI_COMM_WORLD);

                    // compute variance
                    calcGroupQuant(dataMatCopy, val2, valVec2, k, "avg", isField);

                }
                else if (scaleMethod == "sphere")
                {
                    val1 = 0.0;
                    calcGroupQuant(dataMat, val2, valVec2, k, "l2", isField);
                }

                // final calculations for scaling method
                if (isField)
                {
                    if (scaleMethod == "minmax")
                    {
                        for (int i = 0; i < nCells; ++i)
                        {
                            scalingSubVec[k * nCells + i] = valVec1[i];
                            scalingDivVec[k * nCells + i] = valVec2[i] - valVec1[i];
                        }
                    }
                    else if (scaleMethod == "standardize")
                    {
                        for (int i = 0; i < nCells; ++i)
                        {
                            scalingSubVec[k * nCells + i] = valVec1[i];
                            scalingDivVec[k * nCells + i] = sqrt(valVec2[i]);
                        }
                    }
                }
                else
                {
                    if (scaleMethod == "minmax")
                    {
                        scalingSubVec[k] = val1;
                        scalingDivVec[k] = val2 - val1;
                    }
                    else if (scaleMethod == "standardize")
                    {
                        scalingSubVec[k] = val1;
                        scalingDivVec[k] = sqrt(val2);
                    }
                    else if (scaleMethod == "sphere")
                    {
                        scalingSubVec[k] = val1;
                        scalingDivVec[k] = val2;
                    }
                }

            }
            else
            {
                scalingSubVec[k] = 0.0;
                scalingDivVec[k] = scalingInput[k];
            }

            if (!isField)
            {
                cout << "Subtractive scaling factor for " << varName[k] << " is: " << setprecision(numeric_limits<double>::digits10) << scalingSubVec[k] << endl;
                cout << "Divisive scaling factor for " << varName[k] << " is: " << setprecision(numeric_limits<double>::digits10) << scalingDivVec[k] << endl;
            }
            else
            {
                cout << "Scaling field calculation for " << varName[k] << " complete" << endl;
            }

            // distribute to fields in same group instead of recalculating
            for (int g = 0; g < numVars; g++)
            {
                if (g == k)
                {
                    continue;
                }
                else if (scalingInput[g] == scalingInput[k])
                {
                    if (isField)
                    {
                        for (int i = 0; i < nCells; ++i)
                        {
                            scalingSubVec[g * nCells + i] = scalingSubVec[k * nCells + i];
                            scalingDivVec[g * nCells + i] = scalingDivVec[k * nCells + i];
                        }
                    }
                    else
                    {
                        scalingSubVec[g] = scalingSubVec[k];
                        scalingDivVec[g] = scalingDivVec[k];
                    }
                    skipFlags.push_back(k);
                }
            }
        }
        if (scaleMethod == "standardize")
        {
            destroyPMat(dataMatCopy);
        }
    }

    // update global param
    scalingIsField = isField;

    // expand constants to full field for convenience sake
    if (!isField)
    {
        vector<double> subVals(numVars);
        vector<double> divVals(numVars);
        copy(scalingSubVec.begin(), scalingSubVec.end(), subVals.begin());
        copy(scalingDivVec.begin(), scalingDivVec.end(), divVals.begin());
        scalingSubVec.resize(nPoints, 0.0);
        scalingDivVec.resize(nPoints, 0.0);
        for (int k = 0; k < numVars; ++k) 
        {
            for (int i = 0; i < nCells; ++i)
            {
                scalingSubVec[k * nCells + i] = subVals[k];
                scalingDivVec[k * nCells + i] = divVals[k];
            }
        }
    }

    // prevent divisive factors close to zero
    // really only happens when min and max are same, implying field is uniform
    for (int k = 0; k < numVars; ++k) 
    {
        for (int i = 0; i < nCells; ++i)
        {
            if (abs(scalingDivVec[k * nCells + i]) < 1e-15)
            {
                scalingSubVec[k * nCells + i] = 0.0;
                scalingDivVec[k * nCells + i] = 1.0;
            }
        }
    }

    // write scaling fields to file
    cout << "Scaling calculated" << endl;
    if (dataMat->pG->rank == 0)
    {
        // TODO: get SZPLT to output correctly, without silly fileID requirement
        // writeSingle(0, centerVec.data(), "centerProf");
        writeASCIIDoubleVec("scalingSubProf.dat", scalingSubVec);
        writeASCIIDoubleVec("scalingDivProf.dat", scalingDivVec);
    }

    // if also centering, combine into single profile (useful for GEMS)
    if (isCentered)
    {
        cout << "Writing combined subtractive scaling factors" << endl;
        scalingSubVecFull.resize(nPoints);
        for (int i = 0; i < nPoints; ++i)
        {
            scalingSubVecFull[i] = centerVec[i] + scalingSubVec[i];
        }
        if (dataMat->pG->rank == 0)
        {
            writeASCIIDoubleVec("scalingSubProfFull.dat", scalingSubVecFull);
        }
    }

    cout << "Scaling files written" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
}

void tecIO::scaleData(pMat *dataMat)
{
    scaleData(dataMat, false);
}

void tecIO::scaleData(pMat *dataMat, bool unscale)
{

    cout << "Scaling matrix by scaling factor" << endl;
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
                    // (un)scale data
                    if (unscale)
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] *= scalingDivVec[i];
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] += scalingSubVec[i];
                    }
                    else
                    {
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] -= scalingSubVec[i];
                        dataMat->dataD[currentCol * dataMat->myRC[0] + xi + li * dataMat->mb] /= scalingDivVec[i];
                    }
                }
            }
            currentCol++;
        }
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

void tecIO::readDATToVec(std::string filename, std::vector<double> &vec)
{
    if (vec.size() != nPoints)
    {
        cout << "Allocating vector as " << nPoints << " cells" << endl;
        vec.resize(nPoints, 0.0);
        cout << "Vector allocated" << endl;
    }

    // open file, check it exists
    ifstream inFile;
    inFile.open(filename);
    if (!inFile)
    {
        cout << "Failed to open file at " << filename << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    // assumed to be in cell_id order, need hash
    genHash(prefix + to_string(snap0) + suffix);

    // header
    string header;
    inFile >> header;

    // read doubles from file
    int count = 0;
    int varNum;
    double num;
    //vector<double> temp;
    //if (reorder)
    //{
    //    temp.resize(nCells, 0.0);
    //}
    while (inFile >> num)
    {
        
        varNum = count / nCells;
        vec[varNum * nCells + cellID[count % nCells]] = num;
        //if (reorder)
        //{
        //    temp[count % nCells] = num;
        //    if ((count % nCells) == (nCells - 1))
        //    {
        //        for (int i = 0; i < nCells; ++i)
        //        {
        //            vec[varNum * nCells + i] = temp[idx[i]];
        //        }
        //    }
        //}
        count++;
    }
}

void tecIO::readSZPLTToVec(string filename, vector<double> &vec)
{

    if (vec.size() != nPoints)
    {
        cout << "Allocating vector as " << nPoints << " cells" << endl;
        vec.resize(nPoints, 0.0);
        cout << "Vector allocated" << endl;
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    void *fH;
    tecFileReaderOpen(filename.c_str(), &fH);
    int type;
    vector<float> get;
    int ii;
    for (int i = 0; i < numVars; i++)
    {
        ii = getVariableIndex(varName[i], filename);
        tecZoneVarGetType(fH, 1, ii, &type);
        if (type == 1)
        {
            get.resize(nCells);
            tecZoneVarGetFloatValues(fH, 1, ii, 1, nCells, get.data());
            for (int j = 0; j < nCells; j++)
            {
                vec[j + i * nCells] = (double)get[j];
            }
            get.clear();
        }
        else if (type == 2)
        {
            tecZoneVarGetDoubleValues(fH, 1, ii, 1, nCells, &(vec[i * nCells]));
        }
        if (reorder)
        {
            vector<double> temp(nCells, 0.0);
            for (int j = 0; j < nCells; j++)
            {
                temp[j] = vec[i * nCells + j];
            }
            for (int j = 0; j < nCells; j++)
            {
                vec[i * nCells + j] = temp[idx[j]];
            }
            temp.clear();
        }
    }
    tecFileReaderClose(&fH);
    cout << "Vector loaded from: " << filename << endl;
}

void tecIO::vecToCellIDOrder(vector<double> &vecIn, vector<double> &vecOut)
{
    genHash(prefix + to_string(snap0) + suffix);

    //for (int j = 0; j < 
    //vector<double> temp(nCells, 0.0);
    for (int i = 0; i < numVars; ++i)
    {
        for (int j = 0; j < nCells; ++j)
        {
            //temp[j] = vecIn[j * nCells + i];
            vecOut[i * nCells + j] = vecIn[i * nCells + idx[j]];
        }
        //for (int j = 0; j < nCells; ++j)
        //{
        //    vecOut[j * nCells + i] = temp[idx[j]];
        //}
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
