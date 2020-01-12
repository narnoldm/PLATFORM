#include <vector>
#include <string>
#include <iostream>
#include "pMat.hpp"
#include "TECIO.h"

using namespace::std;

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
	void writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat,int nF);
	void writeData(std::string dir, std::string fpref, int wTec, int wBin, int wAscii, pMat *writeMat,int start,int end,int skip);
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
