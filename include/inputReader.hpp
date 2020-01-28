

#ifndef INPUTREADER_H
#define INPUTREADER_H


#include "sysInfo.hpp"
#include "inputInfo.hpp"
#include "operation.hpp"



class inputReader
{
public:
    std::string ifile;
    std::vector<std::string> keys;
    std::vector<int> intParam;
    std::vector<double> doubleParam;
    std::vector<std::string> stringParam;

    sysInfo sys;
    inputInfo inp;
    operationQueue op;
    inputInfo out;

    inputReader(std::string file);
    ~inputReader();
    bool ScanFile();
};


#endif