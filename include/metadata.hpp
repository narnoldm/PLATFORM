
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
    long nPoints, nSets;
    string prefix, suffix;
    vector<string> token;

    meta();
    meta(int t0, int tf, int ts, string &iPrefix, string &iSuffix);
    ~meta();

    virtual void checkSize();
    virtual void checkExists();
    virtual bool readSingle(int fileID, double *point);
    virtual bool writeSingle(int fileID, double *point);
    virtual void miscProcessing(pMat *Mat);
    bool batchWrite(pMat *loadMat);
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

    int numVars, dim, nCells;


    tecIO(int t0, int tf, int ts, string &iPrefix, string &iSuffix, vector<string> &iToken);
    tecIO(tecIO *old, string &dir, string &name);
    ~tecIO();
    virtual void checkSize();
    virtual void checkExists();
    virtual bool readSingle(int fileID, double *point);
    virtual bool writeSingle(int fileID, double *point);
    virtual void miscProcessing(pMat *Mat);

    void addVar(std::string var, string &norm);
    int getVariableIndex(std::string var, std::string file);
    void getDimNodes();
    void genHash();
    void genHash(std::string map);
    void normalize(pMat *dataMat);
    void unNormalize(pMat *dataMat);
    void calcNorm(pMat *dataMat);
    void subAvg(pMat *dataMat);
    void addAvg(pMat *dataMat);
    void calcAvg(pMat *dataMat);
};
