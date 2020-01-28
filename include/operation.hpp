
#ifndef OPERATION_H
#define OPERATION_H

#include<vector>
#include<string>
#include "operation.hpp"
#include "dataID.hpp"
#include "inputInfo.hpp"

/**
 * Operation is defined by one line in input file
 * Operation class also stores all the define operations (gotta be a better way to do this)
 * */
class operation
{
public:

    ///Input names
    std::vector<std::string> input;
    ///Output names
    std::vector<std::string> output;
    ///Input data pointers
    std::vector<dataID *> inputMat;
    ///Output data pointers
    std::vector<dataID *> outputMat;

    /// Operation name
    std::string opName;

    /// Do we know the opName or not
    bool operation_ident;

    operation();
    /// Take token set before equal sign and after it and assign input and output 
    operation(std::string token1, std::string token2);
    ~operation();


    /// Take token and assign input name(s)
    void assignInput(std::string token);
    /// Take token and assign output name(s)
    void assignOutput(std::string token);
    /// Check that names exist in input list and assign associated pointers
    void checkMats(inputInfo &MatList);
    

    /// Loops though inputs to make sure they are defined then assigns output dimensions
    bool inferDim();
    /// Assigns output dimension
    void assignOutputDim();
    /// checks that operation has correct number of inputs and outputs
    void opInfoCheck(const int &,const int &);
    /// Attempts to execute operation
    void execute();
};
/// std out overload
std::ostream &operator<<(std::ostream &os, const operation &op);


/***
 * The operation queue has the set of operations from input file
 * */
class operationQueue
{
public:
    std::vector<operation *> taskQueue;

    operationQueue();
    ~operationQueue();

    ///takes line from parser and adds appropriate operation info
    bool ScanOperations(std::string token);
    /// calls sub operation checkMats to assign pointers
    void assignMats(inputInfo &MatList);
    /// Attempts to infer dimensions of output structures from operation and inputs
    void inferDim();
};

/// std out overload
std::ostream &operator<<(std::ostream &os, const operationQueue &op);
#endif