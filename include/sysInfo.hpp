


#ifndef SYSINFO_H
#define SYSINFO_H


#include "paramID.hpp"


class sysInfo
{
public:
    enum Debug
    {
        IGNORE,
        STOP,
        QUERY
    };
    paramID<long> MEMAVAIL, PC;
    paramID<int> DEBUG;
    sysInfo();
    ~sysInfo();
    bool ScanInfo(string token);
};


std::ostream &operator<<(std::ostream &os, const sysInfo &sID);


#endif