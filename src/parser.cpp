
#include "parser.hpp"

template <typename T>
ostream &operator<<(std::ostream &os, const paramID<T> &pID)
{
        cout << pID.key << " = " << pID.value << " and is of type: " << typeid(pID.value).name() << endl;
        return os;
}
template <class T>
paramID<T>::paramID()
{ 
        key = "";
        value = 0;
        set = false;
};
template <class T>
paramID<T>::paramID(string kIn, T vIn)
{
        key = kIn;
        value = vIn;
        set = true;
};

template <class T>
paramID<T>::paramID(string kIn)
{
        key = kIn;
        value = 0;
        set = false;
};

template <class T>
paramID<T>::~paramID(){};

template <class T>
bool paramID<T>::checkInfo(string line)
{
        int eq = line.find("=");
        if (line.substr(0, eq) != key)
                return 0;
        else
        {
                if (typeid(T) == typeid(int))
                        value = stoi(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(double))
                        value = stod(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(float))
                        value = stof(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(long))
                        value = stol(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(bool))
                        value = to_bool(line.substr(eq + 1));
                else
                {
                        cout << "string to type not defined for " << typeid(T).name() << endl;
                        return 0;
                }
                set = 1;
                return 1;
        }
}
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

dataID::dataID()
{
}

dataID::~dataID()
{
}

void dataID::setName(string n)
{
        name = n;
}

void dataID::setInfo(string d)
{
        token.clear();
        tokenparse(d, "|", token);
        if (token.size() == 1)
        {
                if (d == "inferred")
                {
                        type = inferred;
                        return;
                }
        }
        else //check token 1
        {
                if (token[0] == "input")
                {
                        type = input;
                        return;
                }
                if (token[0] == "defined")
                {
                        type = defined;
                        return;
                }
                if (token[0] == "output")
                {
                        type = output;
                        return;
                }
                cout << token[0] << ": token not recognized" << endl;
                throw(-1);
        }
}
void dataID::setInfo(vector<int> &d)
{
        assert(d.size() > 0);
        dim = d.size();
        dims.resize(d.size());
        for (int i = 0; i < d.size(); i++)
                dims[i] = d[i];
        type = defined;
}
void dataID::switchPmatType(PGrid *newPG)
{
        if (newPG != pMatpoint->pG)
        {
                cout << "in place copy this can be memory explosive" << endl;
                cout << dims[0] << " " << dims[1] << endl;
                pMat *temppMat;
                temppMat = new pMat(dims[0], dims[1], newPG, 0, 0, 0.0);
                temppMat->changeContext(pMatpoint);
                delete pMatpoint;
                pMatpoint = temppMat;
        }
        else
                cout << "already in this format" << endl;
}

ostream &operator<<(std::ostream &os, const dataID &mID)
{
        cout << "Matrix " << mID.name << " is of type:" << mID.type << endl;
        if (mID.type != mID.defined)
        {
                cout << "token : ";
                for (int i = 0; i < mID.token.size(); i++)
                        cout << mID.token[i] << " ";
                cout << endl;
        }
        else
        {
                cout << "dim is: " << mID.dim << " " << endl;
                for (int i = 0; i < mID.dim; i++)
                        cout << "dim[" << i << "]" << mID.dims[i] << " ";
                cout << endl;
                cout << "is Input " << mID.isInput << endl;
                cout << "is Synched Data " << mID.isPA << endl;
                cout << "is P0 Data " << mID.isP0 << endl;
        }
        cout << endl;
        return os;
}

inputInfo::inputInfo()
{
        cout << "setting up mat info" << endl;
}

inputInfo::~inputInfo()
{
        cout << "destructing mat info" << endl;
}

bool inputInfo::ScanInput(string line)
{
        vector<string> split;
        tokenparse(line, "=", split);
        assert(split.size() == 2);
        dataID *point;
        point = new dataID();
        point->setName(split[0]);
        point->setInfo(split[1]);
        matList.push_back(point);
        return 1;
}

bool inputInfo::checkMats()
{
        bool t = true;
        for (int i = 0; i < matList.size(); i++)
        {
                if (matList[i]->type != defined)
                {
                        cout << matList[i]->name << " not defined" << endl;
                        t = false;
                }
        }
        return t;
}

void inputInfo::assignInputInfo()
{
        for (int i = 0; i < matList.size(); i++)
        {
                if (matList[i]->type == 1)
                {
                        cout << (*matList[i]);
                        cout << matList[i]->name << " is an input of type ";
                        cout << matList[i]->token[1] << endl;

                        if (matList[i]->token[1] == "binary")
                        {
                                cout << "input type defined" << endl;
                                pMat temp;
                                matList[i]->dims.resize(2);
                                if (temp.check_bin_size(matList[i]->token[2], matList[i]->dims[0], matList[i]->dims[1]))
                                {
                                        matList[i]->type = defined;
                                        matList[i]->isInput = true;
                                        if (matList[i]->dims[0] == 1 & matList[i]->dims[1] == 1)
                                                matList[i]->dim = scal;
                                        if (matList[i]->dims[0] == 1 & matList[i]->dims[1] == 1)
                                                matList[i]->dim = vec;
                                        else
                                                matList[i]->dim = mat;

                                        cout << (*matList[i]);
                                }
                                else
                                {
                                        cout << "check size failed" << endl;
                                }
                        }
                        if (matList[i]->token[1] == "binaryset")
                        {
                                cout << "input type defined binary set" << endl;
                                matList[i]->datasetInfo = new meta(matList[i]->token);
                                matList[i]->isInput = true;
                                matList[i]->type = defined;
                                matList[i]->dims.resize(2);
                                matList[i]->dims[0] = matList[i]->datasetInfo->nPoints;
                                matList[i]->dims[1] = matList[i]->datasetInfo->nSets;
                                cout << matList[i]->dims[0] << " " << matList[i]->dims[1] << endl;

                                if (matList[i]->dims[0] == 1 & matList[i]->dims[1] == 1)
                                        matList[i]->dim = scal;
                                if (matList[i]->dims[0] == 1 & matList[i]->dims[1] == 1)
                                        matList[i]->dim = vec;
                                else
                                        matList[i]->dim = mat;
                                cout << (*matList[i]);
                        }
                        if (matList[i]->token[1] == "tecplot")
                        {
                                cout << "input type defined tecplot set" << endl;
                                matList[i]->datasetInfo = new tecIO(matList[i]->token);
                                matList[i]->isInput = true;
                                matList[i]->type = defined;
                                matList[i]->dims.resize(2);
                                matList[i]->dims[0] = matList[i]->datasetInfo->nPoints;
                                matList[i]->dims[1] = matList[i]->datasetInfo->nSets;
                                cout << matList[i]->dims[0] << " " << matList[i]->dims[1] << endl;

                                if (matList[i]->dims[0] == 1 & matList[i]->dims[1] == 1)
                                        matList[i]->dim = scal;
                                if (matList[i]->dims[0] == 1 & matList[i]->dims[1] == 1)
                                        matList[i]->dim = vec;
                                else
                                        matList[i]->dim = mat;
                                cout << (*matList[i]);
                        }

                        if (matList[i]->type != defined)
                        {
                                cout << "input not recognized" << endl;
                                throw(-1);
                        }
                }
        }
}

ostream &operator<<(std::ostream &os, const inputInfo &iID)
{
        cout << iID.matList.size() << " data sets defined" << endl;
        for (int i = 0; i < iID.matList.size(); i++)
                cout << (*(iID.matList[i]));
        cout << endl;
        return os;
}

operation::operation(string token1, string token2)
{
        assignInput(token2);
        assignOutput(token1);
}

operation::~operation()
{
}

void operation::assignInput(string token)
{
        input.clear();
        inputMat.clear();
        opName.clear();
        vector<string> split;
        tokenparse(token, "|", split);
        input.resize(split.size() - 1);
        inID.resize(split.size() - 1);
        inputMat.resize(split.size() - 1);
        opName = split[0];
        for (int i = 1; i < split.size(); i++)
                input[i - 1] = split[i];
}

void operation::assignOutput(string token)
{
        output.clear();
        outputMat.clear();
        vector<string> split;
        tokenparse(token, "|", split);
        output.resize(split.size());
        outID.resize(split.size());
        outputMat.resize(split.size());
        for (int i = 0; i < split.size(); i++)
                output[i] = split[i];
}

void operation::checkMats(inputInfo &iInfo)
{
        for (int i = 0; i < inputMat.size(); i++)
        {
                for (int j = 0; j < iInfo.matList.size(); j++)
                {
                        if (iInfo.matList[j]->name == input[i])
                        {
                                inputMat[i] = iInfo.matList[j];
                                inID[i] = j;
                                break;
                        }
                        if (j == iInfo.matList.size() - 1)
                        {
                                cout << input[i] << "name not found" << endl;
                                break;
                        }
                }
        }
        for (int i = 0; i < outputMat.size(); i++)
        {
                for (int j = 0; j < iInfo.matList.size(); j++)
                {
                        if (iInfo.matList[j]->name == output[i])
                        {
                                outputMat[i] = iInfo.matList[j];
                                outID[i] = j;
                                break;
                        }
                        if (j == iInfo.matList.size() - 1)
                        {
                                cout << "name not found" << endl;
                                break;
                        }
                }
        }
}
bool operation::inferDim()
{
        //check all inputs are defined
        for (int i = 0; i < inputMat.size(); i++)
        {
                if (inputMat[i]->type != 3)
                {
                        cout << "input " << i << " " << (*inputMat[i]) << " not defined" << endl;
                        return 0;
                }
        }
        cout << "all inputs defined" << endl;
        assignOutputDim();
        return 1;
}
void operation::assignOutputDim()
{

        if (opName == "svd")
        {
                cout << opName << " recognized" << endl;
                if (inputMat.size() != 1 || outputMat.size() != 3)
                {
                        cout << "SVD should have 3 outputs and 1 input" << endl;
                        assert(inputMat.size() != 1 || outputMat.size() != 3);
                }
                //input contraints SVD requires mb=nb
                inputMat[0]->compPGreq = true;
                //output contraints
                vector<int> temp;
                temp.clear();
                temp.push_back(inputMat[0]->dims[0]);
                temp.push_back(inputMat[0]->dims[1]);
                outputMat[0]->setInfo(temp);
                temp[0] = inputMat[0]->dims[1];
                outputMat[2]->setInfo(temp);
                temp.clear();
                temp.push_back(min(inputMat[0]->dims[0], inputMat[0]->dims[1]));
                outputMat[1]->setInfo(temp);
                outputMat[1]->isPA = true;
                temp.clear();
                return;
        }
        
        else
                cout << opName << " not found" << endl;
}

void operation::execute()
{
        if (opName == "svd")
        {
                cout << "executing SVD" << endl;
                pMat *A, *U, *S, *VT;
                A = inputMat[0]->pMatpoint;
                U = outputMat[0]->pMatpoint;
                S = outputMat[1]->pMatpoint;
                VT = outputMat[2]->pMatpoint;

                A->svd_run(A->N, A->M, 0, 0, U, VT, S->dataD);
                cout << "SVD is destructive: Using or Outputing A is ill advised" << endl;
        }
}
ostream &operator<<(std::ostream &os, const operation &op)
{
        cout << "Operation: " << op.opName << " has " << op.input.size() << " input(s)" << endl;
        for (int i = 0; i < op.input.size(); i++)
                cout << op.input[i] << " " << op.inputMat[i] << " ";
        cout << endl;
        cout << "and " << op.output.size() << " output(s)" << endl;
        for (int i = 0; i < op.output.size(); i++)
                cout << op.output[i] << " " << op.outputMat[i] << " ";
        cout << endl;
}

operationQueue::operationQueue()
{
        cout << "Contructing Operation Queue" << endl;
}
operationQueue::~operationQueue()
{
        cout << "destructing Operation Queue" << endl;
}

bool operationQueue::ScanOperations(string token)
{
        vector<string> split;
        tokenparse(token, "=", split);
        if (split.size() != 2)
                return 0;
        operation *point;
        taskQueue.push_back(point);
        taskQueue[taskQueue.size() - 1] = new operation(split[0], split[1]);
        return 1;
}

void operationQueue::assignMats(inputInfo &MatList)
{
        for (int i = 0; i < taskQueue.size(); i++)
                taskQueue[i]->checkMats(MatList);
}
void operationQueue::inferDim()
{
        for (int i = 0; i < taskQueue.size(); i++)
                taskQueue[i]->inferDim();
}
ostream &operator<<(std::ostream &os, const operationQueue &op)
{
        cout << "Operation Queue has " << op.taskQueue.size() << " operation(s)" << endl;
        for (int i = 0; i < op.taskQueue.size(); i++)
        {
                cout << (*(op.taskQueue[i]));
        }
}

inputReader::inputReader(string file)
{
        cout << "file is :" << file << endl;
        ifile = file;
        ScanFile();
}
inputReader::~inputReader()
{
        ifile.clear();
        keys.clear();
        intParam.clear();
        doubleParam.clear();
        stringParam.clear();
}
bool inputReader::ScanFile()
{
        fstream inFile;
        inFile.open(ifile.c_str(), std::fstream::in);
        string line;
        assert(inFile.good());

        enum blocks
        {
                info,
                input,
                operation,
                output
        };
        int cBlock = info;
        int comment;
        while (!inFile.eof())
        {
                getline(inFile, line);
                comment = line.find("!");
                if (comment != string::npos)
                        line = line.substr(0, comment);
                if (line.empty())
                {
                        continue;
                }
                if (line == "#INFO")
                {
                        cout << line << endl;
                        continue;
                }
                if (line == "#INPUT")
                {
                        cout << line << endl;
                        cBlock = input;
                        continue;
                }
                if (line == "#OPERATIONS")
                {
                        cout << line << endl;
                        cBlock = operation;
                        continue;
                }
                if (line == "#OUTPUT")
                {
                        cout << line << endl;
                        cBlock = output;
                        continue;
                }
                switch (cBlock)
                {
                case info:
                        assert(sys.ScanInfo(line));
                        break;
                case input:
                        assert(inp.ScanInput(line));
                        break;
                case operation:
                        assert(op.ScanOperations(line));
                        break;
                case output:
                        assert(out.ScanInput(line));
                        break;
                }
        }
        cout << sys;
        cout << inp;
        cout << op;
        op.assignMats(inp);
        cout << op;
        cout << out;

        int i = 0;
        while (!inp.checkMats())
        {
                cout << "not all mats defined" << endl;
                inp.assignInputInfo();
                op.inferDim();

                i++;
                if (i > 10)
                {
                        cout << "matricies undefined after 10 iterations: check input file" << endl;
                        throw(-1);
                }
        }
        cout << inp;
        cout << "Yay all orperations and matracies defined" << endl;
}

executioner::executioner()
{
}
executioner::executioner(inputReader *iF)
{
        inpFile = iF;
}
executioner::~executioner()
{
        cout << "destructing Executioner" << endl;
}

void executioner::init()
{
        cout << "initializing matracies" << endl;
        create_matricies();
        cout << "initialization done" << endl;
}
void executioner::exec_all()
{
        cout << "begining execution" << endl;
        for (int i = 0; i < inpFile->op.taskQueue.size(); i++)
        {
                cout << "executing operation " << i << " " << inpFile->op.taskQueue[i]->opName << endl;
                exec(i);
        }
        cout << "execution done" << endl;
}

void executioner::exec(int opID)
{
        inpFile->op.taskQueue[opID]->execute();
}

void executioner::output()
{
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        for (int i = 0; i < inpFile->out.matList.size(); i++)
        {
                if (inpFile->out.matList[i]->token[0] == "output")
                {
                        cout << "output detected" << endl;
                        if (inpFile->out.matList[i]->token[1] == "binary")
                        {
                                cout << "binary detected" << endl;
                                if (rank == 0)
                                        system(("mkdir " + inpFile->out.matList[i]->name).c_str());
                                for (int j = 0; j < inpFile->inp.matList.size(); j++)
                                {
                                        if (inpFile->inp.matList[j]->name == inpFile->out.matList[i]->token[2])
                                        {
                                                cout << inpFile->out.matList[i]->token[2] << " found" << endl;
                                                inpFile->inp.matList[j]->pMatpoint->write_bin((inpFile->out.matList[i]->name + "/" + inpFile->out.matList[i]->name + ".bin").c_str());
                                        }
                                }
                        }
                        if (inpFile->out.matList[i]->token[1] == "binaryset")
                        {
                                cout << "binary set detected" << endl;
                                for (int j = 0; j < inpFile->inp.matList.size(); j++)
                                {
                                        if (inpFile->inp.matList[j]->name == inpFile->out.matList[i]->token[2])
                                        {
                                                cout << inpFile->out.matList[i]->token[2] << " found" << endl;
                                                inpFile->inp.matList[j]->switchPmatType(pGs[1]);
                                                if (inpFile->inp.matList[j]->datasetInfo == NULL)
                                                {
                                                        inpFile->out.matList[i]->datasetInfo = new meta();
                                                        inpFile->out.matList[i]->datasetInfo->snap0 = 1;
                                                        inpFile->out.matList[i]->datasetInfo->snapF = inpFile->inp.matList[j]->pMatpoint->M;
                                                        inpFile->out.matList[i]->datasetInfo->snapSkip = 1;
                                                        inpFile->out.matList[i]->datasetInfo->nPoints = inpFile->inp.matList[j]->pMatpoint->N;
                                                        inpFile->out.matList[i]->datasetInfo->nSets = inpFile->inp.matList[j]->pMatpoint->M;
                                                        inpFile->out.matList[i]->datasetInfo->prefix = inpFile->out.matList[i]->name;
                                                        inpFile->out.matList[i]->datasetInfo->suffix = ".bin";
                                                        inpFile->out.matList[i]->datasetInfo->isInit = true;
                                                }
                                                else
                                                {
                                                        inpFile->out.matList[i]->datasetInfo = inpFile->inp.matList[j]->datasetInfo;
                                                }
                                                inpFile->out.matList[i]->datasetInfo->batchWrite(inpFile->inp.matList[j]->pMatpoint, inpFile->out.matList[i]->name, inpFile->out.matList[i]->name);
                                        }
                                }
                        }
                        if (inpFile->out.matList[i]->token[1] == "tecplot")
                        {
                                cout << "tecplot set detected" << endl;
                                for (int j = 0; j < inpFile->inp.matList.size(); j++)
                                {
                                        if (inpFile->inp.matList[j]->name == inpFile->out.matList[i]->token[2])
                                        {
                                                cout << inpFile->out.matList[i]->token[2] << " found" << endl;
                                                inpFile->inp.matList[j]->switchPmatType(pGs[1]);
                                                tecIO *tempPoint;
                                                if (inpFile->inp.matList[j]->datasetInfo == NULL)
                                                {
                                                        tempPoint = dynamic_cast<tecIO *>(inpFile->out.matList[i]->datasetInfo);
                                                        tempPoint = new tecIO();
                                                        tempPoint->snap0 = 1;
                                                        tempPoint->snapF = inpFile->inp.matList[j]->pMatpoint->M;
                                                        tempPoint->snapSkip = 1;
                                                        tempPoint->nSets = inpFile->inp.matList[j]->pMatpoint->M;
                                                        tempPoint->prefix = inpFile->out.matList[i]->name;
                                                        tempPoint->suffix = ".szplt";
                                                        tempPoint->isInit = true;
                                                        tempPoint->meshFile = inpFile->out.matList[i]->token[3];
                                                        tempPoint->fixedMesh = true;
                                                        tempPoint->getDimNodes();
                                                        for (int k = 4; k < inpFile->out.matList[i]->token.size(); k += 2)
                                                        {
                                                                tempPoint->addVarO(inpFile->out.matList[i]->token[k], inpFile->out.matList[i]->token[k + 1]);
                                                        }
                                                        tempPoint->nPoints = tempPoint->nCells * tempPoint->numVars;
                                                }
                                                else
                                                {
                                                        tempPoint = dynamic_cast<tecIO *>(inpFile->inp.matList[j]->datasetInfo);
                                                }

                                                tempPoint->batchWrite(inpFile->inp.matList[j]->pMatpoint, inpFile->out.matList[i]->name, inpFile->out.matList[i]->name);
                                                delete tempPoint;
                                        }
                                }
                        }
                }
        }
}
void executioner::clear()
{
        for (int i = 0; i < inpFile->inp.matList.size(); i++)
        {
                delete inpFile->inp.matList[i]->pMatpoint;
                cout << inpFile->inp.matList[i]->name << endl;
        }
        for (int i = 0; i < pGs.size(); i++)
        {
                delete pGs[i];
        }
}
void executioner::create_matricies()
{
        int rank, size;
        PGrid *temp;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        temp = new PGrid(rank, size, 0);
        pGs.push_back(temp);
        temp = new PGrid(rank, size, 1);
        pGs.push_back(temp);
        pMat *pointMat;
        //pMats.clear();
        for (int i = 0; i < inpFile->inp.matList.size(); i++)
        {
                cout << "allocating matrix " << *(inpFile->inp.matList[i]) << endl;
                if (inpFile->inp.matList[i]->isPA)
                {
                        cout << "Creating Synched data" << endl;
                        pointMat = new pMat(inpFile->inp.matList[i]->dims[0], pGs[1]->size, pGs[1], 0, 0, 0.0);
                }
                else if (inpFile->inp.matList[i]->isInput)
                {
                        cout << "Creating Loading matrix" << endl;
                        pointMat = new pMat(inpFile->inp.matList[i]->dims[0], inpFile->inp.matList[i]->dims[1], pGs[1], 0, 0, 0.0);
                }
                else
                        pointMat = new pMat(inpFile->inp.matList[i]->dims[0], inpFile->inp.matList[i]->dims[1], pGs[0], 0, 0, 0.0);
                //point header to pMat
                inpFile->inp.matList[i]->pMatpoint = pointMat;

                //input
                if (inpFile->inp.matList[i]->isInput)
                {
                        cout << "loading " << *(inpFile->inp.matList[i]) << endl;
                        if (inpFile->inp.matList[i]->token[1] == "binary")
                        {
                                pointMat->read_bin(inpFile->inp.matList[i]->token[2]);
                        }
                        else if (inpFile->inp.matList[i]->token[1] == "binaryset")
                        {
                                inpFile->inp.matList[i]->datasetInfo->batchRead(pointMat);
                        }
                        else if (inpFile->inp.matList[i]->token[1] == "tecplot")
                        {
                                tecIO *tempPoint = dynamic_cast<tecIO *>(inpFile->inp.matList[i]->datasetInfo);
                                tempPoint->batchRead(pointMat);
                        }
                        else
                        {
                                cout << "load type unrecognized" << endl;
                                throw(-1);
                        }
                }

                if (inpFile->inp.matList[i]->compPGreq)
                {
                        cout << "need to copy to even process grid" << endl;
                        inpFile->inp.matList[i]->switchPmatType(pGs[0]);
                }
                //pMats.push_back(pointMat);
        }
}

ostream &operator<<(std::ostream &os, const executioner &e)
{
}

bool to_bool(std::string str)
{
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        std::istringstream is(str);
        bool b;
        is >> std::boolalpha >> b;
        return b;
}

void tokenparse(string &input, string sep, vector<string> &tokens)
{
        //count |
        int loc, loc2;
        tokens.clear();
        loc = input.find(sep);
        if (loc == string::npos)
        {
                tokens.push_back(input);
                return;
        }
        else
        {
                tokens.push_back(input.substr(0, loc));
                while (loc != string::npos)
                {
                        loc2 = loc;
                        loc = input.find(sep, loc2 + 1);
                        if (loc == string::npos)
                        {
                                tokens.push_back(input.substr(loc2 + 1));
                                break;
                        }
                        else
                        {
                                tokens.push_back(input.substr(loc2 + 1, loc - loc2 - 1));
                        }
                }
        }
}