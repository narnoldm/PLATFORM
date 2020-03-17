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
    assert(snapF > snap0);
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
            fid = fopen((prefix + to_string(i) + suffix).c_str(), "rb");
            assert(fid != NULL);
            fclose(fid);
        }
    }
}

bool meta::readSingle(int fileID, double *point)
{
    cout << "meta read" << endl;
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
    assert(loadMat->mb == nPoints);
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
    assert(system(NULL)); //check if system commands work
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (!rank)
        system(("mkdir " + dir).c_str());
    assert(loadMat->mb == nPoints);
    int iP = 0, fileIndex, localC = 0;
    if (isInit)
        fileIndex = snap0;
    else
    {
        fileIndex = 1;
        nSets = loadMat->N;
    }
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

            cout << "proc " << iP << " is writing file " << fileIndex << endl;
            writeSingle(fileIndex, loadMat->dataD.data() + nPoints * localC, dir + "/" + fpref);
            localC++;
        }
    }
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
    assert(snapF > snap0);
    nSets = 0;
    for (int i = snap0; i <= snapF; i = i + snapSkip)
        nSets++;
    checkSize();
    checkExists();
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
    cout << (prefix + std::to_string(fileID) + suffix) << endl;
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
        tecFileReaderOpen((prefix + std::to_string(fileID) + suffix).c_str(), &infH);
    assert(infH != NULL);
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
        tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, jMax, &point[(i - dim) * jMax]);
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
    cout << "metadata has registered " << numVars << " variable(s)" << endl;
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
            delete vName;
            vName = NULL;
        }
        delete vName;
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
    //synch nodes and dimension
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

void tecIO::genHash()
{
    int var_index = 0;
        var_index = getVariableIndex("cell_id");
        hash.resize(nCells, 0);
        cellID.resize(nCells, 0);
        void *fH;
        if (printRank)
        {
                int hashType;
                cellID.resize(nCells);
                tecFileReaderOpen((prefix + std::to_string(snap0) + suffix).c_str(), &fH);
                tecZoneVarGetType(fH, 1, var_index, &hashType);
                if (hashType == 3)
                {
                        tecZoneVarGetInt32Values(fH, 1, var_index, 1, nCells, cellID.data());
                }
                if (hashType == 2)
                {
                        std::vector<float> hashTemp(nCells, 0.0);
                        tecZoneVarGetFloatValues(fH, 1, var_index, 1, nCells, hashTemp.data());
                        for (int i = 0; i < hash.size(); i++)
                                hash[i] = (int)(hashTemp[i]);
                        hashTemp.clear();
                }
                tecFileReaderClose(&fH);
                for (int i = 0; i < nCells; i++)
                {
                        cellID[i]--;
                        hash[cellID[i]] = i;
                }
                printf("hash table built\nDistributing to procs\n");
        }
        MPI_Bcast(hash.data(), nCells, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(cellID.data(), nCells, MPI_INT, 0, MPI_COMM_WORLD);

}

void tecIO::genHash(string map)
{
}

void tecIO::normalize(pMat *dataMat)
{
    cout<<"Normalizing Matrix by Norm Factor"<<endl;
        int numFiles = dataMat->nelements / nPoints;
        for (int i = 0; i < numFiles; i++)
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

void tecIO::unNormalize(pMat *dataMat)
{
    cout<<"UnNormalizing Matrix by Norm Factor"<<endl;
        int numFiles = dataMat->nelements / nPoints;
        for (int i = 0; i < numFiles; i++)
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

void tecIO::calcNorm(pMat *dataMat)
{
     int numFiles = dataMat->nelements / nPoints;
        std::vector<double> mag(nCells * numFiles, 0.0);
        //double *mag = new double[nCells*numFiles];
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
                        
                        cout<<varName[i]<<" is in normalization group "<< normFactor[i]<<endl;
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
                                        cout<<"Setting norm factor "<<j<<" to "<<maxmag<<endl;
                                }
                        }
                }
                else
                {
                        cout<<varName[i]<<" is using specified normalization of "<< normFactor[i];
                }
        }
        mag.clear();
        for (int i = 0; i < numVars; i++)
        {
                MPI_Allreduce(MPI_IN_PLACE, normFactor.data() + i, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
                cout<<"Synched norm factor for "<<varName[i]<<" is : "<< normFactor[i];
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
                printASCIIVecP0("norm.data",normVec.data(),normVec.size());
                normVec.clear();
                //delete [] normVec;
        }



}

void tecIO::subAvg(pMat *dataMat)
{
    if (average.size() == 0)
    {
        cout<<"Average isn't setup call calcAvg first"<<endl;
        throw(-1);
    }
    
    cout<<"Subtracting Average"<<endl;
    int numFiles = dataMat->nelements / nPoints;
    for (int i = 0; i < numFiles; i++)
    {
        for (int j = 0; j < nPoints; j++)
            dataMat->dataD[i * nPoints + j] -= average[j];
    }
}

void tecIO::calcAvg(pMat *dataMat)
{
    if (average.size() != nPoints)
    {
        cout << "Allocating Average as " << nPoints << " cells" << endl;
        average.resize(nPoints, 0.0);
        cout << "Average Allocated" << endl;
    }
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
        writeSingle(snap0, average.data(), "average");
}
