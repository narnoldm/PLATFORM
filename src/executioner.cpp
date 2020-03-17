#include "executioner.hpp"


using namespace::std;


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
                                                inpFile->inp.matList[j]->switchPmatType(1);
                                                if (inpFile->inp.matList[j]->datasetInfo == NULL)
                                                {
                                                        inpFile->out.matList[i]->datasetInfo = new meta();
                                                        inpFile->out.matList[i]->datasetInfo->snap0 = 1;
                                                        inpFile->out.matList[i]->datasetInfo->snapF = inpFile->inp.matList[j]->pMatpoint->N;
                                                        inpFile->out.matList[i]->datasetInfo->snapSkip = 1;
                                                        inpFile->out.matList[i]->datasetInfo->nPoints = inpFile->inp.matList[j]->pMatpoint->M;
                                                        inpFile->out.matList[i]->datasetInfo->nSets = inpFile->inp.matList[j]->pMatpoint->N;
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
                                                inpFile->inp.matList[j]->switchPmatType(1);
                                                tecIO *tempPoint;
                                                if (inpFile->inp.matList[j]->datasetInfo == NULL)
                                                {
                                                        tempPoint = dynamic_cast<tecIO *>(inpFile->out.matList[i]->datasetInfo);
                                                        tempPoint = new tecIO();
                                                        tempPoint->snap0 = 1;
                                                        tempPoint->snapF = inpFile->inp.matList[j]->pMatpoint->N;
                                                        tempPoint->snapSkip = 1;
                                                        tempPoint->nSets = inpFile->inp.matList[j]->pMatpoint->N;
                                                        tempPoint->prefix = inpFile->out.matList[i]->name;
                                                        tempPoint->suffix = ".szplt";
                                                        tempPoint->isInit = true;
                                                        tempPoint->meshFile = inpFile->out.matList[i]->token[3];
                                                        tempPoint->fixedMesh = true;
                                                        tempPoint->getDimNodes();
                                                        for (int k = 4; k < inpFile->out.matList[i]->token.size(); k += 2)
                                                        {
                                                                tempPoint->addVar(inpFile->out.matList[i]->token[k], inpFile->out.matList[i]->token[k + 1]);
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
        pMat *pointMat;
        //pMats.clear();
        for (int i = 0; i < inpFile->inp.matList.size(); i++)
        {
                cout << "allocating matrix " << *(inpFile->inp.matList[i]) << endl;
                if (inpFile->inp.matList[i]->isPA)
                {
                        cout << "Creating Synched data" << endl;
                        pointMat = new pMat(inpFile->inp.matList[i]->dims[0]*pGs[0]->prow, pGs[0]->pcol, pGs[0], 0, 3, 0.0);
                }
                else if (inpFile->inp.matList[i]->isInput)
                {
                        cout << "Creating Loading matrix" << endl;
                        pointMat = new pMat(inpFile->inp.matList[i]->dims[0], inpFile->inp.matList[i]->dims[1], pGs[0], 0, 1, 0.0);
                }
                else
                        pointMat = new pMat(inpFile->inp.matList[i]->dims[0], inpFile->inp.matList[i]->dims[1], pGs[0], 0, 0, 0.0);
                //point header to pMat
                inpFile->inp.matList[i]->pMatpoint = pointMat;

                cout<<(*pointMat);

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
                        inpFile->inp.matList[i]->switchPmatType(0);
                }
                //pMats.push_back(pointMat);
        }
}

ostream &operator<<(std::ostream &os, const executioner &e)
{
}
