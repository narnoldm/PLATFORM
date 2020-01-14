#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <type_traits>
#include <typeinfo>
#include <assert.h>
#include <algorithm>
#include <cctype>

#include "pMat.hpp"

using namespace ::std;

template <class T>
class paramID
{
public:
    string key;
    bool set;
    T value;

    paramID();
    paramID(string kIn, T vIn);
    paramID(string kIn);
    ~paramID();

    bool checkInfo(string line);
};

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
    vector<string> token;
    string name;
    vector<int> dims;

    dataID();
    ~dataID();
    void setName(string n);
    void setInfo(string d);
    void setInfo(vector<int> &d);
};

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
    vector<dataID *> matList;

    inputInfo();
    ~inputInfo();

    bool ScanInput(string token);
    bool checkMats();
    void assignInputInfo();
};

class operation
{
public:

    vector<string> input;
    vector<string> output;
    vector<int> inID;
    vector<int> outID;
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
};

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
    bool ScanFile();
};


class executioner
{
    public: 
    inputReader * inpFile;
    //vector<dataTool*> dTs;
    vector<PGrid*> pGs;
    vector<pMat*> pMats;


    executioner();
    executioner(inputReader * iF);
    ~executioner();

    void init();
    void exec_all();
    void exec(int opID);
    void clear();

    void create_matricies();
};



bool to_bool(std::string str);

void tokenparse(string &input, string sep, vector<string> &tokens);