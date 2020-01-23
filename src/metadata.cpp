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
        for (int i = snap0; i < snapF; i = i + snapSkip)
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
    //checkExists();
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
}
bool tecIO::readSingle(int fileID, double *point)
{
}
bool tecIO::writeSingle(int fileID, double *point)
{
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

dataTool::dataTool(int r, std::string p, std::string s, int t0, int ts, int tf)
{
    prefix = p;
    suffix = s;
    snap0 = t0;
    snapF = tf;
    numVars = 0;
    snapSkip = ts;
    rank = r;
    printRank = (rank == 0);
    nCells = 0;
    N = 0;
    M = 0;
}

dataTool::~dataTool()
{
    if (printRank)
        printf("Deconstructing data Tool");
}

void dataTool::genHash()
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
void dataTool::normalize(pMat *dataMat)
{
    if (printRank)
        printf("Normalizing Matrix by Norm Factor\n");
    int numFiles = dataMat->nelements / N;
    for (int i = 0; i < numFiles; i++)
    {
        for (int j = 0; j < numVars; j++)
        {
            for (int k = 0; k < nCells; k++)
            {
                dataMat->dataD[i * N + j * nCells + k] = dataMat->dataD[i * N + j * nCells + k] / normFactor[j];
            }
        }
    }
}
void dataTool::unnormalize(pMat *dataMat)
{
    if (printRank)
        printf("Normalizing Matrix by Norm Factor\n");
    int numFiles = dataMat->nelements / N;
    for (int i = 0; i < numFiles; i++)
    {
        for (int j = 0; j < numVars; j++)
        {
            for (int k = 0; k < nCells; k++)
            {
                dataMat->dataD[i * N + j * nCells + k] = dataMat->dataD[i * N + j * nCells + k] * normFactor[j];
            }
        }
    }
}
void dataTool::calcNorm(pMat *dataMat)
{
    int numFiles = dataMat->nelements / N;
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
            if (printRank)
                printf("%s is in normalization group %f\n", varName[i].c_str(), -normFactor[i]);
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
                    if (printRank)
                        printf("Setting norm factor %d to %f\n", j, maxmag);
                }
            }
        }
        else
        {
            if (printRank)
                printf("%s is using specified normalization of %f\n", varName[i].c_str(), normFactor[i]);
        }
    }
    mag.clear();
    for (int i = 0; i < numVars; i++)
    {
        MPI_Allreduce(MPI_IN_PLACE, normFactor.data() + i, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if (printRank)
        {
            printf("Synched norm factor for %s is : %f\n", varName[i].c_str(), normFactor[i]);
        }
    }
    if (printRank)
    {
        std::vector<double> normVec(numVars * nCells);
        for (int i = 0; i < numVars; i++)
        {
            for (int j = 0; j < nCells; j++)
            {
                normVec[i * nCells + j] = normFactor[i];
            }
        }
        writeTec("norm", rank + 1, normVec.data());
        writeAscii("norm", rank + 1, normVec.data());
        normVec.clear();
        //delete [] normVec;
    }
}
/*void dataTool::calcError(pMat *&ROMMat, pMat *&projFOM, pMat *&error, std::string filename)
{
        if (printRank)
                printf("caclulating ROM Error\n");
        int numFiles = ROMMat->nelements / N;
        if (error->M != M)
        {
                std::cout << "dimension mismatch" << std::endl;
                throw(-1);
        }
        for (int i = 0; i < numFiles; i++)
        {
                error->dataD[i] = cblas_dnrm2(N, &(ROMMat->dataD[i * N]), 1) / cblas_dnrm2(N, &(projFOM->dataD[i * N]), 1);
        }
        error->write_bin(filename);
}*/

void dataTool::subAvg(pMat *&dataMat)
{
    if (average.size() == 0)
    {
        if (printRank)
            printf("Average isn't setup call calcAvg first\n");
        throw(-1);
    }
    if (printRank)
        printf("Subtracting Average\n");
    int numFiles = dataMat->nelements / N;
    for (int i = 0; i < numFiles; i++)
    {
        for (int j = 0; j < N; j++)
            dataMat->dataD[i * N + j] = dataMat->dataD[i * N + j] - average[j];
    }
}

void dataTool::addAvg(pMat *&dataMat)
{
    if (average.size() == 0)
    {
        if (printRank)
            printf("Average isn't setup call calcAvg first\n");
        throw(-1);
    }
    if (printRank)
        printf("Subtracting Average\n");
    int numFiles = dataMat->nelements / N;
    for (int i = 0; i < numFiles; i++)
    {
        for (int j = 0; j < N; j++)
            dataMat->dataD[i * N + j] = dataMat->dataD[i * N + j] + average[j];
    }
}

void dataTool::calcAvg(pMat *dataMat, int wTec, int wAscii)
{
    if (average.size() != N)
    {
        if (printRank)
            printf("Allocating Average as %d cells\n", N);
        average.resize(N, 0.0);
        if (printRank)
            printf("Average Allocated\n");
    }
    int numFiles = dataMat->nelements / N;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < numFiles; j++)
        {
            average[i] += dataMat->dataD[j * N + i];
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, average.data(), N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    for (int i = 0; i < N; i++)
        average[i] = average[i] / M;
    if (printRank && (dim != 1))
    {
        printf("Average calculated\n");
        if (wTec == 1)
            writeTec("average", rank + 1, average.data());
        if (wAscii == 1)
            writeAscii("average", rank + 1, average.data());
    }
}
void dataTool::readAvg(std::string filename)
{
    if (average.size() != N)
    {
        if (printRank)
            printf("Allocating Average as %d cells\n", N);
        average.resize(N, 0.0);
        if (printRank)
            printf("Average Allocated\n");
    }

    if (rank == 0)
    {
        void *fH;
        tecFileReaderOpen(filename.c_str(), &fH);
        int type;
        std::vector<float> get;
        for (int i = 0; i < numVars; i++)
        {
            tecZoneVarGetType(fH, 1, i + dim + 1, &type);
            if (type == 1)
            {
                //get= new float[nCells];
                get.resize(nCells);
                tecZoneVarGetFloatValues(fH, 1, i + dim + 1, 1, nCells, get.data());
                for (int j = 0; j < nCells; j++)
                {
                    average[j + i * nCells] = (double)get[j];
                }
                //delete [] get;
                get.clear();
            }
            else if (type == 2)
            {
                tecZoneVarGetDoubleValues(fH, 1, i + dim + 1, 1, nCells, &(average[i * nCells]));
            }
        }
        tecFileReaderClose(&fH);
        std::cout << "average loaded from :" << filename << std::endl;
    }

    MPI_Bcast(average.data(), N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (rank == 0)
        std::cout << "average Broad-casted" << std::endl;
}
void dataTool::writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat)
{
    writeData(dir, fpref, wTec, wBin, wAscii, writeMat, 0, writeMat->M, 1);
}
void dataTool::writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat, int nF)
{
    writeData(dir, fpref, wTec, wBin, wAscii, writeMat, 0, nF, 1);
}
void dataTool::writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat, int start, int end, int skip)
{
    //checks
    MPI_Barrier(MPI_COMM_WORLD);
    if (printRank)
        system(("mkdir " + dir).c_str());

    MPI_Barrier(MPI_COMM_WORLD);
    if (writeMat->pG->prow != 1)
    {
        printf("blocking must be rectangular for benign I/O\n");
        throw(-1);
    }
    if (writeMat->nelements % N != 0)
    {
        printf("distribute array and tecplot mismatch\n");
        throw(-1);
    }
    int numFiles = writeMat->nelements / N;
    int iP = 0;
    int fileIndex = 0;
    int localC = 0;
    for (int i = start; i < end; i = i + skip)
    {
        iP = (int)(i / writeMat->mb);
        while (iP > (writeMat->pG->size - 1))
        {
            iP = iP - writeMat->pG->size;
        }
        if (writeMat->pG->rank == iP)
        {
            fileIndex = i + 1;
            if (printRank)
                printf("proc %d is writing file %d\n", iP, fileIndex);
            if (wTec == 1)
                writeTec(dir + "/" + fpref, fileIndex, writeMat->dataD.data() + N * localC);
            if (wAscii == 1)
                writeAscii(dir + "/" + fpref, fileIndex, writeMat->dataD.data() + N * localC);
            if (wBin == 1)
            {
                write_single_bin(dir + "/" + fpref, fileIndex, writeMat->dataD.data() + N * localC);
            }
            localC++;
        }
    }
}

void dataTool::writeAscii(std::string fpref, int fileIndex, double *data)
{
    FILE *fid;
    if (hash.size() == 0)
    {
        printf("hash table not setup exiting\n");
        throw(-1);
    }
    if ((fid = fopen((fpref + std::to_string(fileIndex) + ".dat").c_str(), "w")) == NULL)
    {
        printf("error with file open\n");
    }
    fprintf(fid, "variables= %s\n", (fpref + std::to_string(fileIndex)).c_str());
    for (int i = 0; i < numVars; i++)
    {
        for (int j = 0; j < nCells; j++)
            fprintf(fid, "%16.16E\n", data[i * nCells + hash[j]]);
    }
}

int dataTool::write_single_bin(std::string fpref, int fileIndex, double *data)
{
    FILE *fid;
    fid = fopen((fpref + std::to_string(fileIndex) + ".bin").c_str(), "wb");
    int ONE = 1;
    fwrite(&N, sizeof(int), 1, fid);
    fwrite(&ONE, sizeof(int), 1, fid);
    for (int i = 0; i < numVars; i++)
    {
        for (int j = 0; j < nCells; j++)
            fwrite(&(data[i * nCells + hash[j]]), sizeof(double), 1, fid);
    }
    fclose(fid);
    return 1;
}
int dataTool::write_single_pbin(std::string fpref, int fileIndex, double *data)
{
    FILE *fid;
    fid = fopen((fpref + std::to_string(fileIndex) + ".pbin").c_str(), "wb");
    int ONE = 1;
    fwrite(&N, sizeof(int), 1, fid);
    fwrite(&ONE, sizeof(int), 1, fid);
    for (int i = 0; i < numVars; i++)
    {
        for (int j = 0; j < nCells; j++)
            fwrite(&(data[i * nCells + hash[j]]), sizeof(double), 1, fid);
    }
    //fwrite(data,sizeof(double),N,fid);
    fclose(fid);
    return 1;
}

void dataTool::writeTec(std::string fpref, int fileIndex, double *data)
{
    void *infH = NULL;
    void *outfH = NULL;
    tecFileReaderOpen((prefix + std::to_string(snap0) + suffix).c_str(), &infH);
    int zoneType;
    long iMax, jMax, kMax;
    tecZoneGetType(infH, 1, &zoneType);
    tecZoneGetIJK(infH, 1, &iMax, &jMax, &kMax);
    if ((zoneType != 5) && (zoneType != 3))
    {
        printf("Zone is weird\n");
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

    /*int *varTypes = new int[dim+numVars];
                int *valueLoc = new int[dim+numVars];
                int *passive = new int[dim+numVars];
                int *shareVar = new int[dim+numVars];*/
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
    for (int i = 0; i < dim; i++)
    {
        if (varTypes[i] == 1) //float
        {
            float *nDat;
            if (i == 0)
                nDat = new float[iMax];
            tecZoneVarGetFloatValues(infH, 1, i + 1, 1, iMax, &nDat[0]);
            tecZoneVarWriteFloatValues(outfH, 1, i + 1, 0, iMax, &nDat[0]);
            if (i == (dim - 1))
                delete[] nDat;
        }
        else if (varTypes[i] == 2)
        {
            double *nDat;
            if (i == 0)
                nDat = new double[iMax];
            tecZoneVarGetDoubleValues(infH, 1, i + 1, 1, iMax, &nDat[0]);
            tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, iMax, &nDat[0]);
            if (i == (dim - 1))
                delete[] nDat;
        }
    }
    for (int i = dim; i < (dim + numVars); i++)
    {
        tecZoneVarWriteDoubleValues(outfH, 1, i + 1, 0, jMax, &data[(i - dim) * jMax]);
    }
    long numValues;
    tecZoneNodeMapGetNumValues(infH, 1, jMax, &numValues);
    int *nodeMap = new int[numValues];
    tecZoneNodeMapGet(infH, 1, 1, jMax, &nodeMap[0]);
    tecZoneNodeMapWrite32(outfH, 1, 0, 1, numValues, &nodeMap[0]);
    delete[] nodeMap;

    tecFileReaderClose(&infH);
    tecFileWriterClose(&outfH);
    varTypes.clear();
    valueLoc.clear();
    passive.clear();
    shareVar.clear();
    /*delete [] varTypes;
                delete [] valueLoc;
                delete [] passive;
                delete [] shareVar;*/
}
void dataTool::loadDataTec(pMat *loadMat)
{
    //checks
    if (loadMat->pG->prow != 1)
    {
        printf("blocking must be rectangular for benign read in\n");
        throw(-1);
    }
    if (loadMat->nelements % N != 0)
    {
        printf("distribute array and tecplot mismatch\n");
        throw(-1);
    }
    int numFiles = loadMat->nelements / N;
    int iP = 0;
    int fileIndex = 0;
    int localC = 0;
    for (int i = 0; i < M; i++)
    {
        iP = (int)(i / loadMat->mb);
        while (iP > (loadMat->pG->size - 1))
        {
            iP = iP - loadMat->pG->size;
        }
        if (loadMat->pG->rank == iP)
        {
            fileIndex = snap0 + i * snapSkip;
            if (printRank)
                printf("proc %d is reading file %d\n", iP, fileIndex);
            load(fileIndex, loadMat->dataD.data() + N * localC);
            localC++;
        }
    }
}

void dataTool::load(int fileIndex, double *data)
{
    void *fH;
    tecFileReaderOpen((prefix + std::to_string(fileIndex) + suffix).c_str(), &fH);
    int type;
    std::vector<float> get;
    for (int i = 0; i < numVars; i++)
    {
        tecZoneVarGetType(fH, 1, varIndex[i], &type);
        if (type == 1)
        {
            //get= new float[nCells];
            get.resize(nCells);
            tecZoneVarGetFloatValues(fH, 1, varIndex[i], 1, nCells, get.data());
            for (int j = 0; j < nCells; j++)
            {
                data[j + i * nCells] = (double)get[j];
            }
            //delete [] get;
            get.clear();
        }
        else if (type == 2)
        {
            tecZoneVarGetDoubleValues(fH, 1, varIndex[i], 1, nCells, &(data[i * nCells]));
        }
    }
    tecFileReaderClose(&fH);
}

void dataTool::loadFile(const char *filename, void *data, int numVariables, int *variableIndex)
{
    void *fH;
    tecFileReaderOpen(filename, &fH);
    int type = 0;
    std::vector<float> get;
    for (int i = 0; i < numVariables; i++)
    {
        printf("%d\n", variableIndex[i]);
        tecZoneVarGetType(fH, 1, variableIndex[i], &type);
        printf("TYPE: %d\n", type);
        if (type == 1)
        {
            get.resize(nCells);
            tecZoneVarGetFloatValues(fH, 1, variableIndex[i], 1, nCells, get.data());
            double *hdata = (double *)data;
            for (int j = 0; j < nCells; j++)
            {
                hdata[j + i * nCells] = (double)get[j];
            }
            get.clear();
        }
        else if (type == 2)
        {
            double *hdata = (double *)data;
            tecZoneVarGetDoubleValues(fH, 1, variableIndex[i], 1, nCells, &(hdata[i * nCells]));
        }
        else if (type == 3)
        {
            int *hdata = (int *)data;
            tecZoneVarGetInt32Values(fH, 1, variableIndex[i], 1, nCells, &(hdata[i * nCells]));
        }
    }
    tecFileReaderClose(&fH);
}

void dataTool::loadDataBin(pMat *loadMat)
{
    //checks
    if (loadMat->pG->prow != 1)
    {
        printf("blocking must be rectangular for benign read in\n");
        throw(-1);
    }
    int numFiles = loadMat->nelements / N;
    int iP = 0;
    int fileIndex = 0;
    int localC = 0;
    for (int i = 0; i < M; i++)
    {
        iP = (int)(i / loadMat->mb);
        while (iP > (loadMat->pG->size - 1))
        {
            iP = iP - loadMat->pG->size;
        }
        if (loadMat->pG->rank == iP)
        {
            fileIndex = snap0 + i * snapSkip;
            if (printRank)
                printf("proc %d is reading file %d\n", iP, fileIndex);
            read_single_bin(fileIndex, loadMat->dataD.data() + N * localC);
            localC++;
        }
    }
}

int dataTool::read_single_bin(int index, double *data)
{
    FILE *fid;
    char *id = new char[6];
    sprintf(id, "%06d", index);
    //printf("checking %s\n",id);
    fid = fopen((prefix + std::to_string(index) + suffix).c_str(), "rb");
    int Mfile, Nfile;
    fread(&Nfile, sizeof(int), 1, fid);
    fread(&Mfile, sizeof(int), 1, fid);
    if (Mfile != 1)
    {
        printf("file has more than 1 snapshot !!!\n");
        return -1;
    }
    if (Nfile != N)
    {
        printf("file has incorrect number of points");
        return -1;
    }
    fread(data, sizeof(double), N, fid);
    fclose(fid);
    return 1;
}

void dataTool::prepTool()
{
    if (numVars == 0)
    {
        if (printRank)
            printf("No tec vars selected assuming bin read\n");
        numVars = 1;
    }
    N = nCells * numVars;
    M = 0;
    for (int i = snap0; i <= snapF; i = i + snapSkip)
    {
        M++;
    }
    if (printRank)
    {
        printf("Data Tool preped to load N = %d = %d vars * %d cells , M = %d snapshots\n", N, numVars, nCells, M);
    }
}
void dataTool::addVar(std::string var, double norm)
{
    varName.push_back(var);
    varIndex.push_back(getVariableIndex(var));
    normFactor.push_back(norm);
    if (varName.size() != varIndex.size())
    {
        if (printRank)
        {
            printf("Error in var arrays\n");
        }
        throw(-1);
    }
    numVars = varName.size();
    if (printRank)
        printf("dataTool now has %d variables\n", numVars);
}

int dataTool::getVariableIndex(std::string var)
{
    void *fH;
    int fileVars;
    char *fileName = NULL;
    int tecIndex = 0;
    if (printRank)
    {
        tecFileReaderOpen((prefix + std::to_string(snap0) + suffix).c_str(), &fH);
        tecDataSetGetNumVars(fH, &fileVars);
        for (int i = 1; i <= fileVars; i++)
        {
            tecVarGetName(fH, i, &fileName);
            if (var == fileName)
            {
                printf("%s found in index %d\n", fileName, i);
                tecIndex = i;
                break;
            }
            fileName = NULL;
        }
        tecFileReaderClose(&fH);
        if (tecIndex == 0)
        {
            printf("Var not found\n");
            throw(-1);
        }
    }
    MPI_Bcast(&tecIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return (tecIndex);
}
int dataTool::getVariableIndex(std::string var, std::string file)
{
    void *fH;
    int fileVars;
    char *fileName = NULL;
    int tecIndex = 0;
    if (printRank)
    {
        tecFileReaderOpen(file.c_str(), &fH);
        tecDataSetGetNumVars(fH, &fileVars);
        for (int i = 1; i <= fileVars; i++)
        {
            tecVarGetName(fH, i, &fileName);
            if (var == fileName)
            {
                printf("%s found in index %d\n", fileName, i);
                tecIndex = i;
                break;
            }
            fileName = NULL;
        }
        tecFileReaderClose(&fH);
        if (tecIndex == 0)
        {
            printf("Var not found\n");
            throw(-1);
        }
    }
    MPI_Bcast(&tecIndex, 1, MPI_INT, 0, MPI_COMM_WORLD);
    return (tecIndex);
}

void dataTool::getPointsBin()
{
    FILE *fid;
    int Mfile, Nfile;
    dim = 1;
    //char *id=new char[5];
    if (printRank)
    {
        //printf("%06d\n",snap0);
        //sprintf(id,"%06d",snap0);
        //printf("%s\n",id);
        printf("%s\n", (prefix + std::to_string(snap0) + suffix).c_str());
        fid = fopen((prefix + std::to_string(snap0) + suffix).c_str(), "rb");
        fread(&Nfile, sizeof(int), 1, fid);
        fread(&Mfile, sizeof(int), 1, fid);
        nCells = Nfile;
        fclose(fid);
    }
    MPI_Bcast(&nCells, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dim, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

void dataTool::getDimNodes()
{
    void *fH;
    long iMax = 0, jMax = 0, kMax = 0;
    if (printRank)
    {
        long check;
        int id = 0;
        printf("%s\n", (prefix + std::to_string(snap0) + suffix).c_str());
        tecFileReaderOpen((prefix + std::to_string(snap0) + suffix).c_str(), &fH);
        tecZoneGetIJK(fH, 1, &iMax, &jMax, &kMax);
        printf("iMax: %d, jMax %d, kMax %d \n", iMax, jMax, kMax);
        check = iMax;
        while (check == iMax)
        {
            id++;
            tecZoneVarGetNumValues(fH, 1, id, &check);
            printf("var %d has %d values\n", id, check);
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