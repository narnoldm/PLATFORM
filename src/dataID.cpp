
#include "dataID.hpp"


using namespace::std;

dataID::dataID()
{
}

dataID::~dataID()
{
}

void dataID::setName(string &n)
{
        name = n;
}

void dataID::setInfo(string &d)
{
        token.clear();
        tokenparse(d, "|", token);
        if (token.size() == 1)
        {
                if (token[0] == "inferred")
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
        if(!(d.size()>0) )
                throw(-1);
        dim = d.size();
        dims.resize(d.size());
        for (int i = 0; i < d.size(); i++)
                dims[i] = d[i];
        type = defined;
}
void dataID::switchPmatType(int newblock)
{
        cout<<newblock<<" "<<pMatpoint->block<<endl;
        if (newblock != pMatpoint->block)
        {
                cout << "in place copy this can be memory explosive" << endl;
                cout << dims[0] << " " << dims[1] << endl;
                pMat *temppMat;
                temppMat = new pMat(dims[0], dims[1], pMatpoint->pG, 0, newblock, 0.0);
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