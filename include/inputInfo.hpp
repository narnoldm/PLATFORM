
#ifndef INPUTINFO_H
#define INPUTINFO_H

#include "dataID.hpp"

/**
 * Input info is the parser for elements in the INPUT block of the input file
 * It will attempt to assign dataID's to each of the defined names
 * */
class inputInfo
{
public:
    enum typeenum
    {
        inferred,
        input,
        output,
        defined
    };
    enum dimenum
    {
        scal,
        vec,
        mat,
    };

    ///Pointer set to each of the dataID objects
    std::vector<dataID *> matList;

    inputInfo();
    ~inputInfo();

    /// Parser for individual input line
    bool ScanInput(std::string token);

    ///  Checks to see if the matracies all have defined sizes
    bool checkMats();
    ///  If identified as input will call input reader to determine size
    void assignInputInfo();
};

/// std out overload
std::ostream &operator<<(std::ostream &os, const inputInfo &iID);

#endif