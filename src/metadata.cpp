#include "metadata.hpp"

meta::meta()
{
}
meta::~meta()
{
}

meta::meta(int t0, int tf, int ts, string &iPrefix, string &iSuffix)
{
    snap0 = t0;
    snapF = tf;
    snapSkip = ts;
    prefix = iPrefix;
    suffix = iSuffix;
    assert(snapF > snap0);
    nSets = 0;
    for (int i = snap0; i <= snapF; i = i + snapSkip)
        nSets++;
    cout << prefix + to_string(snap0) + suffix << endl;
    checkSize();
    checkExists();
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
            fid = fopen((prefix + to_string(i) + suffix).c_str(), "rb");
            assert(fid != NULL);
            fclose(fid);
        }
    }
}

bool meta::readSingle(int fileID, double *point)
{
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
    assert(loadMat->pG->prow == 1);
    int iP = 0;
    int fileIndex = snap0;
    int localC = 0;
    for (int i = 0; i < nSets; i++)
    {
        iP = (int)(i / loadMat->mb);
        while (iP > (loadMat->pG->size - 1))
        {
            iP = iP - loadMat->pG->size;
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
}

bool meta::writeSingle(int fileID, double *point)
{
    FILE *fid;
    fid = fopen((prefix + to_string(fileID) + suffix).c_str(), "rb");

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
    assert(loadMat->pG->prow == 1);
    int iP = 0;
    int fileIndex = snap0;
    int localC = 0;
    for (int i = 0; i < nSets; i++)
    {
        iP = (int)(i / loadMat->mb);
        while (iP > (loadMat->pG->size - 1))
        {
            iP = iP - loadMat->pG->size;
        }
        if (loadMat->pG->rank == iP)
        {
            fileIndex = snap0 + i * snapSkip;

            cout << "proc " << iP << " is writing file " << fileIndex << endl;
            writeSingle(fileIndex, loadMat->dataD.data() + nPoints * localC);
            localC++;
        }
    }
}

void meta::miscProcessing(pMat *Mat)
{
    cout << "no additional processing for binary" << endl;
}

tecIO::tecIO(int t0, int tf, int ts, string &iPrefix, string &iSuffix, vector<string> &iToken)
{
    snap0 = t0;
    snapF = tf;
    snapSkip = ts;
    prefix = iPrefix;
    suffix = iSuffix;
    token = iToken;
    assert(snapF > snap0);
    nSets = 0;
    for (int i = snap0; i <= snapF; i = i + snapSkip)
        nSets++;
    cout << prefix + to_string(snap0) + suffix << endl;
    checkSize();
    checkExists();
}

void tecIO::checkSize()
{
    getDimNodes();
    for (int i = 0; i < token.size(); i = i + 2)
    {
        addVar(token[i], token[i + 1]);
    }
    nPoints = nCells * numVars;
}

void tecIO::checkExists()
{
    void *fH = NULL;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
    {
        for (int i = snap0; i <= snapF; i = i + snapSkip)
        {
            tecFileReaderOpen((prefix + std::to_string(i) + suffix).c_str(), &fH);
            assert(fH != NULL);
            tecFileReaderClose(&fH);
        }
    }
}
bool tecIO::readSingle(int fileID, double *point)
{
    void *fH=NULL;
    tecFileReaderOpen((prefix + std::to_string(fileIndex) + suffix).c_str(), &fH);
    int type;
    assert(fH!=NULL);
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
                data[j + i * nCells] = (double)get[j];
            }
            get.clear();
        }
        else if (type == 2)
        {
            tecZoneVarGetDoubleValues(fH, 1, varIndex[i], 1, nCells, &(data[i * nCells]));
        }
    }
    tecFileReaderClose(&fH);
    
}
bool tecIO::writeSingle(int fileID, double *point)
{
    void *infH = NULL;
    void *outfH = NULL;
    tecFileReaderOpen((prefix + std::to_string(snap0) + suffix).c_str(), &infH);
    assert(infH!=NULL);
    int zoneType;
    long iMax, jMax, kMax;
    tecZoneGetType(infH, 1, &zoneType);
    tecZoneGetIJK(infH, 1, &iMax, &jMax, &kMax);
    if ((zoneType != 5) && (zoneType != 3))
    {
        printf("Zone is weird/Not supported\n");
        throw(-1);
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
    tecFileWriterOpen((fpref + std::to_string(fileIndex) + suffix).c_str(), "Code out", varstr.c_str(), 1, 0, 1, NULL, &outfH);
    assert(outfH!=NULL)
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
    tecZoneCreateFE(outfH, std::to_string(fileIndex).c_str(), zoneType, iMax, jMax, &varTypes[0], &shareVar[0], &valueLoc[0], &passive[0], shareCon, 0, 0, &outZone);
    tecZoneSetUnsteadyOptions(outfH, outZone, fileIndex, (int)fileIndex);

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
            tecZoneVarGetDoubleValues(infH, 1, i + 1, 1, iMax, &nDat[0]);
            tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, iMax, &nDat[0]);
        }
    }
    nfDat.clear();
    ndDat.clear();
    for (int i = dim; i < (dim + numVars); i++)
    {
        tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, jMax, &data[(i - dim) * jMax]);
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
void tecIO::miscProcessing(pMat *Mat)
{
}

void tecIO::addVar(string var, string &norm)
{
    varName.push_back(var);
    varIndex.push_back(getVariableIndex(var, prefix + to_string(snap0) + suffix));
    normID.push_back(norm);
    double temp;
    temp = stod(norm);
    normFactor.push_back(temp);
    assert(varName.size() == varIndex.size());
    numVars = varName.size();
    cout << "metadata registers " << numVars << " variables" << endl;
}

int tecIO::getVariableIndex(string var, string file)
{
    void *fH;
    int fileVars;
    char *vName = NULL;
    int tecIndex = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
    {
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
            vName = NULL;
        }
        vName = NULL;
        tecFileReaderClose(&fH);
        if (tecIndex == 0)
        {
            cout << "Var not found :" << var << endl;
            throw(-1);
        }
    }
    MPI_Bcast(&tecIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return (tecIndex);
}

void tecIO::getDimNodes()
{
    void *fH = NULL;
    long iMax = 0, jMax = 0, kMax = 0;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
    {
        long check;
        int id = 0;
        cout << (prefix + std::to_string(snap0) + suffix) << endl;
        tecFileReaderOpen((prefix + std::to_string(snap0) + suffix).c_str(), &fH);
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
    //synch nodes and dimension
    MPI_Bcast(&nCells, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void tecIO::genHash()
{
}

void tecIO::genHash(string map)
{
}

void tecIO::normalize(pMat *dataMat)
{
}

void tecIO::unNormalize(pMat *dataMat)
{
}

void tecIO::calcNorm(pMat *dataMat)
{
}

void tecIO::subAvg(pMat *dataMat)
{
}

void tecIO::calcAvg(pMat *dataMat)
{
}
