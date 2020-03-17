
#ifndef DATAID_H
#define DATAID_H


#include<vector>
#include<string>

#include "metadata.hpp"
#include "misc_func.hpp"

/**
 *  The dataID object is the wrapper for the loaded data.
 *  It contains information about the current known information
 *  This includes if its size is known as well as if it is fully defined
 *  It also keeps track of any operation requirment for the given data
 **/

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
    /// Keeps track of if matrix dimensions are known
    int type, dim;

    /// Tags for if matrix needs to be in input format, synched format, or single processor format
    bool isInput=false,isPA=false,isP0=false;
     
    ///Tags that even blocking is required for operation
    bool compPGreq=false;
    ///Tags that this is an IO matrix
    bool IOPGreq=false;

    ///Pointer to pMat(ScaLAPACK) data structure
    pMat *pMatpoint=NULL;
    ///Pointer to meta data input if input or output
    meta *datasetInfo=NULL;
    ///Token from input block for additional processing
    std::vector<std::string> token;
    ///Structure name
    std::string name;
    ///Dimension of Structure 
    std::vector<int> dims;


    dataID();
    ~dataID();
    /// Sets internal name
    void setName(std::string &n);
    /// Sets intenal dimension sets
    void setInfo(std::string &d);
    /// Sets dimensions
    void setInfo(std::vector<int> &d);
    /// For pMat structures performs processer grid swap needed for IO or computational constraints
    void switchPmatType(int newtype);
};
/// stdout overload
std::ostream &operator<<(std::ostream &os, const dataID &mID);



#endif