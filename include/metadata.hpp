
#include "TECIO.h"
#include <string>
#include <assert.h>
#include <vector>
#include "pMat.hpp"

using namespace ::std;

class meta
{
public:
    int snap0, snapF, snapSkip;
    bool isInit = false;
    long nPoints, nSets;
    string prefix, suffix;
    vector<string> token;

    meta();
    meta(vector<string> &iToken);
    meta(int t0, int tf, int ts, string &iPrefix, string &iSuffix);
    virtual ~meta();
    virtual void init(int t0, int tf, int ts, string &iPrefix, string &iSuffix);

    virtual void checkSize();
    virtual void checkExists();
    virtual bool readSingle(int fileID, double *point);
    virtual bool writeSingle(int fileID, double *point, string fpref);
    virtual void miscProcessing(pMat *Mat);
    bool batchWrite(pMat *loadMat);
    bool batchWrite(pMat *loadMat, string dir, string fpref);
    bool batchRead(pMat *loadMat);
};


class tecIO : public meta
{
public:
    //child specific

    vector<string> varName;
    vector<int> varIndex;
    vector<string> normID;
    vector<double> normFactor;
    vector<int> hash;
    vector<int> cellID;
    vector<double> average;
    string meshFile;
    bool fixedMesh = false;

    int numVars, dim, nCells;

    tecIO();
    tecIO(vector<string> &iToken);
    tecIO(int t0, int tf, int ts, string &iPrefix, string &iSuffix);
    virtual void init(int t0, int tf, int ts, string &iPrefix, string &iSuffix);

    ~tecIO();
    virtual void checkSize();
    virtual void checkExists();
    virtual bool readSingle(int fileID, double *point);
    virtual bool writeSingle(int fileID, double *point, string fpref);
    virtual void miscProcessing(pMat *Mat);

    bool writeSingleFull(int fileID, double *point, string fpref, string meshfile);
    void addVar(std::string var, string &norm);
    int getVariableIndex(std::string var, std::string file);
    void getDimNodes();
    void checkMeshDim(string filename);
    void genHash();
    void genHash(std::string map);
    void normalize(pMat *dataMat);
    void unNormalize(pMat *dataMat);
    void calcNorm(pMat *dataMat);
    void subAvg(pMat *dataMat);
    void addAvg(pMat *dataMat);
    void calcAvg(pMat *dataMat);
};
