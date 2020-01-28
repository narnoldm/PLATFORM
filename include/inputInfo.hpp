
#ifndef INPUTINFO_H
#define INPUTINFO_H


#include "dataID.hpp"

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
    std::vector<dataID *> matList;

    inputInfo();
    ~inputInfo();

    bool ScanInput(std::string token);
    bool checkMats();
    void assignInputInfo();
};


std::ostream &operator<<(std::ostream &os, const inputInfo &iID);


#endif