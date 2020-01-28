#include "inputInfo.hpp"


using namespace::std;


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