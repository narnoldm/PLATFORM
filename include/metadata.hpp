//H**********************************************************************
// FILENAME :        metadata.hpp    
//
// DESCRIPTION :
//       Header file for metadata file information headers
//       Contains Meta class and associated methods
//
// Univeristy of Michigan Ann Arbor 
// Computational AeroSciences Laboratory
// AUTHOR :    Nicholas Arnold-Medabalimi        
// DATE :      6/17/2022
//
//H*/

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




///
/// Meta contains the information about how the data files should
/// be organized in the parallel matrix. Overloads of the metaclass
/// are used to support extra file types.
///
class meta
{
public:
    int snap0, snapF, snapSkip; ///< The indexing of the filenames snap0 through snapF with a skip of snapSkip
    bool isInit = false; ///< Flag to indicate if the meta data has been initialized
    long nPoints, nSets; ///< The number of points and number of files (nSets)
    std::string prefix, suffix; ///< The prefix and suffix of the filenames
    std::vector<std::string> token; ///< The input token used to initalize meta

    meta(); ///< Default constructor
    meta(std::vector<std::string> &iToken);///< Constructor that takes in the input token
    meta(int t0, int tf, int ts, std::string &iPrefix, std::string &iSuffix);///< Constructor of decoupled token
    virtual ~meta();///< Destructor
    virtual void init(int t0, int tf, int ts, std::string &iPrefix, std::string &iSuffix);///< Initialize meta

    virtual void checkSize(); ///< Check the size of a single file
    virtual void checkExists();///< Check if the files exist
    virtual bool readSingle(int fileID, double *point);///< Read a single file
    virtual bool writeSingle(int fileID, double *point, std::string fpref);///< Write a single file
    virtual bool writeSingle(int fileID, double *point, std::string fpref, int points);///< Write a single file if points is diiferent from initalized meta
    virtual void miscProcessing(pMat *Mat);///< Perform misc processing on the data (not used in parent meta)
    bool batchWrite(pMat *loadMat); ///< Write the data to the full set of files using defaults
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref, int nModes); ///< Write a subset of the files
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref, int mStart, int mEnd, int mSkip); ///< Write a subset of the files
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip);///< Write a subset of the files
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref); ///< Write the data to the full set of files using specified names
    bool batchWrite(pMat *loadMat, std::string dir, std::string fpref, int mStart, int mEnd, int mSkip, int fStart, int fSkip, int dim);
    bool batchRead(pMat *loadMat); ///< Read the data from the full set of files using defaults
    bool batchRead(pMat *loadMat, int ii); ///< Read the data from the full set of files using defaults
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
    std::vector<double> centerVec;

    // data scaling
    // (x - scalingSubVec) / scalingDivVec
    bool isScaled = false;
    bool scalingIsField;
    std::vector<double> scalingInput;
    std::vector<double> scalingSubVec, scalingDivVec, scalingSubVecFull;

    std::vector<int> idx;
    std::vector<int> cellID;

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
    virtual void miscProcessing(pMat *Mat);

    // I/O
    void vecToCellIDOrder(std::vector<double> &vecIn, std::vector<double> &vecOut);
    void readDATToVec(std::string filename, std::vector<double> &vec);
    void readSZPLTToVec(std::string filename, std::vector<double> &vec);
    virtual bool readSingle(int fileID, double *point);
    virtual bool writeSingle(int fileID, double *point, std::string fpref);
    virtual bool writeSingle(int fileID, double *point, std::string fpref, int points);
    bool writeSingleFile(std::string filename, std::vector<std::string> &fvars, double *point, std::string meshfile);

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
