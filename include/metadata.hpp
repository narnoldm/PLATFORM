


#ifndef METADATA_H
#define METADATA_H

#include "TECIO.h"
#include <string>
#include <assert.h>
#include <vector>
#include "pMat.hpp"
#include "misc_func.hpp"

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
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref);
    bool batchRead(pMat *loadMat);
};


class tecIO : public meta
{
public:
    //child specific

    std::vector<std::string> varName;
    std::vector<int> varIndex;
    std::vector<std::string> normID;
    std::vector<double> normFactor;
    std::vector<int> hash;
    std::vector<int> cellID;
    std::vector<double> average;
    std::string meshFile;
    bool outBin=false;
    bool fixedMesh = false;
    bool GEMSbin=false;

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

    bool writeSingleFull(int fileID, double *point, std::string fpref, std::string meshfile);
    void addVar(std::string var, std::string &norm);
    int getVariableIndex(std::string var, std::string file);
    void getDimNodes();
    void checkMeshDim(std::string filename);
    void genHash(std::string );
    void normalize(pMat *dataMat);
    void unNormalize(pMat *dataMat);
    void calcNorm(pMat *dataMat);
    void subAvg(pMat *dataMat);
    void addAvg(pMat *dataMat);
    void calcAvg(pMat *dataMat);
    void readAvg(std::string filename);

    //misc
    void activateGEMSbin(std::string);
};


#endif
