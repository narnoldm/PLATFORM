
#ifndef DATAID_H
#define DATAID_H


#include<vector>
#include<string>

#include "metadata.hpp"
#include "misc_func.hpp"



class dataID
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
    int type, dim;
    bool isInput=0,isPA=0,isP0=0;
     
     
    bool compPGreq=false;
    bool IOPGreq=false;


    pMat *pMatpoint=NULL;
    meta *datasetInfo=NULL;
    std::vector<std::string> token;
    std::string name;
    std::vector<int> dims;

    


    dataID();
    ~dataID();
    void setName(std::string n);
    void setInfo(std::string d);
    void setInfo(std::vector<int> &d);
    void switchPmatType(PGrid *newPG);
};

std::ostream &operator<<(std::ostream &os, const dataID &mID);



#endif