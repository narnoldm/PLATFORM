
#ifndef OPERATION_H
#define OPERATION_H

#include<vector>
#include<string>
#include "operation.hpp"
#include "dataID.hpp"
#include "inputInfo.hpp"

using namespace::std;


class operation
{
public:

    vector<string> input;
    vector<string> output;
    vector<dataID *> inputMat;
    vector<dataID *> outputMat;



    string opName;
    bool operation_ident;

    operation();
    operation(string token1, string token2);
    ~operation();

    void assignInput(string token);
    void assignOutput(string token);
    void checkMats(inputInfo &MatList);
    void checkOperation(string token);
    bool inferDim();
    void assignOutputDim();
    void opInfoCheck(int,int);
    void execute();
};
std::ostream &operator<<(std::ostream &os, const operation &op);
class operationQueue
{
public:
    vector<operation *> taskQueue;

    operationQueue();
    ~operationQueue();

    bool ScanOperations(string token);
    void assignMats(inputInfo &MatList);
    void inferDim();
};
std::ostream &operator<<(std::ostream &os, const operationQueue &op);
#endif