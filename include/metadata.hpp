

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
    virtual void readSingle(int fileID, double *point);
    virtual void readSingleLowMem(int fileID, pMat* dataMat, int colIdx);
    virtual void writeSingle(int fileID, double *point, std::string fpref);
    virtual void writeSingle(int fileID, double *point, std::string fpref, int points);
    void batchWrite(pMat *loadMat);
    void batchWrite(pMat *loadMat, std::string dir, std::string fpref, int nModes);
    void batchWrite(pMat *loadMat, std::string dir, std::string fpref, int mStart, int mEnd, int mSkip);
    void batchWrite(pMat *loadMat, std::string dir, std::string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip);
    void batchWrite(pMat *loadMat, std::string dir, std::string fpref);
    void batchWrite(pMat *loadMat, std::string dir, std::string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip, int dim);
    void batchRead(pMat *loadMat);
    void batchRead(pMat *loadMat, int ii);
};

class tecIO : public meta
{
public:
    //child specific

    std::vector<std::string> varName;
    std::vector<int> varIndex;
    std::vector<std::string> normID;

    // data centering
    // x - centerVec
    bool isCentered = false;
    bool centeringIsField;
    std::vector<double> centeringInput;
    pMat* centerVec = NULL;

    // data scaling
    // (x - scalingSubVec) / scalingDivVec
    bool isScaled = false;
    bool scalingIsField;
    std::vector<double> scalingInput;
    pMat* scalingSubVec = NULL;
    pMat* scalingDivVec = NULL;
    pMat* scalingSubVecFull = NULL;

    std::vector<int> idx;

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

    // I/O
    void vecToCellIDOrder(std::vector<double> &vecIn, std::vector<double> &vecOut);
    void readSZPLTToPMat(std::string filename, pMat* loadMat);
    void readDATToVec(std::string filename, std::vector<double> &vec);
    void readSZPLTToVec(std::string filename, std::vector<double> &vec);
    virtual void readSingle(int fileID, double *point);
    virtual void readSingle(std::string filename, double *point);
    virtual void readSingleLowMem(int fileID, pMat* dataMat, int colIdx);
    virtual void readSingleLowMem(std::string filename, pMat* dataMat, int colIdx);
    void writeSingle(int fileID, double *point, std::string fpref);
    void writeSingle(int fileID, double *point, std::string fpref, int points);
    void writeSingleFile(std::string filename, std::vector<std::string> &fvars, double *point, std::string meshfile);

    void addVar(std::string var, std::string &norm);
    int getVariableIndex(std::string var, std::string file);
    void getDimNodes();
    void checkMeshDim(std::string filename);
    void genHash(std::string);

    // feature scaling routines
    void calcGroupQuant(pMat *dataMat, double &outVal, std::vector<double> &outVec, int varIdx, std::string methodName, bool aggCells);
    void calcCentering(pMat *dataMat, std::string centerMethod);
    void calcCentering(pMat *dataMat, std::string centerMethod, bool isField);
    void calcCentering(pMat *dataMat, std::string centerMethod, bool isField, bool writeToDisk);
    void calcScaling(pMat *dataMat, std::string scaleMethod);
    void calcScaling(pMat *dataMat, std::string scaleMethod, bool isField);
    void calcScaling(pMat *dataMat, std::string scaleMethod, bool isField, bool writeToDisk);
    void scaleData(pMat *dataMat);
    void scaleData(pMat *dataMat, bool unscale);
    void centerData(pMat *dataMat);
    void centerData(pMat *dataMat, bool uncenter);

    //misc
    void activateGEMSbin(std::string);
    void activateReorder(std::string);
};

int compareMeta(meta* meta1, meta* meta2);

#endif
