
#include "sysInfo.hpp"

sysInfo::sysInfo()
{
        cout << "constructing system info" << endl;
        MEMAVAIL.key = "MEMMAX";
        PC.key = "PROCSCHECK";
        DEBUG.key = "DEBUG";
}
sysInfo::~sysInfo()
{
        cout << "deconstructing system info" << endl;
}

bool sysInfo::ScanInfo(string token)
{
        MEMAVAIL.checkInfo(token);
        PC.checkInfo(token);
        DEBUG.checkInfo(token);
        return true;
}

ostream &operator<<(std::ostream &os, const sysInfo &sID)
{
        cout << sID.MEMAVAIL << sID.PC << sID.DEBUG;
        return os;
}
