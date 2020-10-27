

#ifndef INPUTREADER_H
#define INPUTREADER_H

#include "sysInfo.hpp"
#include "inputInfo.hpp"
#include "operation.hpp"

/**
 * The input reader class is responsible for reading in parsing the input file into 
 * the varius sub input file objects
 * These include:
 * System information (memory/proc limits) 
 * Input Info (data to be read in and its associated size)
 * Operations (List of operations to be done)
 * Output Info (How resulting data should be output)
 */
class inputReader
{

public:
    /// input file name
    std::string ifile;
    /// Sytem Information
    sysInfo sys;
    /// Data Structure Information
    inputInfo inp;
    /// Operation sets
    operationQueue op;
    /// Output Information
    inputInfo out;

    /// Contstructor will generate the subset
    inputReader(std::string file);
    /// Destructor
    ~inputReader();
    /// Explicit function that scans inputs and populates subset
    bool ScanFile();
};

#endif