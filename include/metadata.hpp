
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

class tecIO : meta
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

class dataTool
{
public:
    int snap0, snapF, snapSkip, numVars, M, dim, rank;
    long nCells, N;
    bool printRank;
    std::string prefix, suffix;
    std::vector<std::string> varName;
    std::vector<int> varIndex;
    std::vector<double> normFactor;

    std::vector<int> hash;
    std::vector<int> cellID;
    std::vector<double> average;

public:
    dataTool(int r, std::string p, std::string s, int t0, int ts, int tf);
    ~dataTool();
    void genHash();
    void genHash(std::string map);
    void normalize(pMat *dataMat);
    void unnormalize(pMat *dataMat);
    void calcNorm(pMat *dataMat);
    //void calcError(pMat *&ROMMat,pMat *&projFOM,pMat *&error,std::string filename);
    void subAvg(pMat *&dataMat);
    void addAvg(pMat *&dataMat);
    void calcAvg(pMat *dataMat, int wTec, int wAscii);
    void readAvg(std::string filename);
    void writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat);
    void writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat, int nF);
    void writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat, int start, int end, int skip);
    void writeAscii(std::string fpref, int fileIndex, double *data);
    int write_single_bin(std::string fpref, int fileIndex, double *data);
    int write_single_pbin(std::string fpref, int fileIndex, double *data);
    void writeTec(std::string fpref, int fileIndex, double *data);
    void loadDataTec(pMat *loadMat);
    void load(int fileIndex, double *data);
    void loadFile(const char *filename, void *data, int numVariables, int *variableIndex);
    void loadDataBin(pMat *loadMat);
    int read_single_bin(int index, double *data);
    void prepTool();
    void addVar(std::string var, double norm);
    int getVariableIndex(std::string var);
    int getVariableIndex(std::string var, std::string file);
    void getPointsBin();
    void getDimNodes();
};
