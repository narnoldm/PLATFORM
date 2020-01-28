
#include "operation.hpp"


using namespace::std;

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
                                break;
                        }
                        if (j == iInfo.matList.size() - 1)
                        {
                                cout << input[i] << "name not found" << endl;
                                throw(-1);
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
                                break;
                        }
                        if (j == iInfo.matList.size() - 1)
                        {
                                cout << "name not found" << endl;
                                throw(-1);
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
                opInfoCheck(1,3);
                //input contraints SVD requires mb=nb
                inputMat[0]->compPGreq = true;
                //output contraints
                int minMN= min(inputMat[0]->dims[0], inputMat[0]->dims[1]);
                vector<int> temp;
                temp.clear();
                temp.push_back(inputMat[0]->dims[0]);
                temp.push_back(minMN);
                outputMat[0]->setInfo(temp);
                temp[0] = minMN;
                temp[1] = inputMat[0]->dims[1];
                outputMat[2]->setInfo(temp);


                temp.clear();
                temp.push_back(minMN);
                outputMat[1]->setInfo(temp);
                outputMat[1]->isPA = true;
                temp.clear();
                return;
        }
        else if (opName == "M*M")
        {
                cout << opName << " recognized" << endl;
                opInfoCheck(2,1);
                //input contraints SVD requires mb=nb
                inputMat[0]->compPGreq = true;
                inputMat[1]->compPGreq = true;
                //output contraints
                vector<int> temp;
                temp.clear();
                temp.push_back(inputMat[0]->dims[0]);
                temp.push_back(inputMat[1]->dims[1]);
                outputMat[0]->setInfo(temp);
                return;
        }
        else if (opName == "Transpose")
        {
                cout<<opName<<"recognized"<<endl;
                opInfoCheck(1,1);

        }
        else
        {
                cout << opName << " not found" << endl;
                throw(-1);
        }
}
void operation::execute()
{
        if (opName == "svd")
        {
                cout << "executing"<< opName<< endl;
                pMat *A, *U, *S, *VT;
                A = inputMat[0]->pMatpoint;
                U = outputMat[0]->pMatpoint;
                S = outputMat[1]->pMatpoint;
                VT = outputMat[2]->pMatpoint;

                A->svd_run(A->N, A->M, 0, 0, U, VT, S->dataD);
                cout << "SVD is destructive: Using or Outputing A is ill advised" << endl;
        }
        else if(opName =="M*M")
        {
                cout << "executing "<<opName << endl;
                pMat *C,*A,*B; // C= A*B;
                A=inputMat[0]->pMatpoint;
                B=inputMat[1]->pMatpoint;
                C=outputMat[0]->pMatpoint;
                assert(inputMat[0]->dims[1]==inputMat[1]->dims[0]);
                C->matrix_Product('N','N',outputMat[0]->dims[0],outputMat[0]->dims[1],inputMat[0]->dims[1],inputMat[0]->pMatpoint,0,0,inputMat[1]->pMatpoint,0,0,1.0,0.0,0,0);
        }

}

void operation::opInfoCheck(const int &in,const int &out)
{
        if (inputMat.size() != in || outputMat.size() != out)
        {
                cout << opName <<" should have "<<out<<" output(s) and "<<in<<" input(s)" << endl;
                assert(inputMat.size() != in || outputMat.size() != out);
        }
        return;
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
        point=new operation(split[0], split[1]);
        taskQueue.push_back(point);
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