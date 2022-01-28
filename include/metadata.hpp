

#ifndef METADATA_H
#define METADATA_H

#include "TECIO.h"
#include <string>
#include <numeric>
#include <assert.h>
#include <vector>
#include "pMat.hpp"
#include "misc_func.hpp"
#include "misc_func_mpi.hpp"

class meta
{
public:
    int snap0, snapF, snapSkip;
    bool isInit = false;
    long nPoints, nSets;
    std::string prefix, suffix;
    std::vector<std::string> token;

    meta();
    meta(std::vector<std::string> &iToken);
    meta(int t0, int tf, int ts, std::string &iPrefix, std::string &iSuffix);
    virtual ~meta();
    virtual void init(int t0, int tf, int ts, std::string &iPrefix, std::string &iSuffix);

    virtual void checkSize();
    virtual void checkExists();
    virtual bool readSingle(int fileID, double *point);
    virtual bool writeSingle(int fileID, double *point, std::string fpref);
    virtual void miscProcessing(pMat *Mat);
    bool batchWrite(pMat *loadMat);
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref, int);
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref, int, int, int);
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref);
    bool batchRead(pMat *loadMat);
    bool batchRead(pMat *loadMat, int ii);
};

class tecIO : public meta
{
public:
    //child specific

    std::vector<std::string> varName;
    std::vector<int> varIndex;
    std::vector<std::string> normID;
    std::vector<double> normFactor;
    std::vector<int> idx;
    std::vector<int> cellID;
    std::vector<double> average;
    std::string meshFile;
    bool outBin = false;
    bool fixedMesh = false;
    bool GEMSbin = false;
    bool reorder = false;

    int numVars, dim;
    long nCells;

    tecIO();
    tecIO(std::vector<std::string> &iToken);
    tecIO(int t0, int tf, int ts, std::string &iPrefix, std::string &iSuffix);
    virtual void init(int t0, int tf, int ts, std::string &iPrefix, std::string &iSuffix);

    ~tecIO();
    virtual void checkSize();
    virtual void checkExists();
    virtual bool readSingle(int fileID, double *point);
    virtual bool writeSingle(int fileID, double *point, std::string fpref);
    virtual void miscProcessing(pMat *Mat);

    bool writeSingleFile(std::string filename, std::vector<std::string> &fvars, double *point, std::string meshfile);
    void addVar(std::string var, std::string &norm);
    int getVariableIndex(std::string var, std::string file);
    void getDimNodes();
    void checkMeshDim(std::string filename);
    void genHash(std::string);
    void normalize(pMat *dataMat);
    void unNormalize(pMat *dataMat);
    void calcNorm(pMat *dataMat);
    void subAvg(pMat *dataMat);
    void addAvg(pMat *dataMat);
    void calcAvg(pMat *dataMat);
    void readAvg(std::string filename);

    //misc
    void activateGEMSbin(std::string);
    void activateReorder(std::string);
};

int compareMeta(meta* meta1, meta* meta2);

#endif
