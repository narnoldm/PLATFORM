
#include "misc_func.hpp"


using namespace::std;

void printASCIIVecP0(std::string fname, double *Mat, int N)
{
        FILE *fid;
        fid = fopen(fname.c_str(), "w");
        int i = 0;
        for (i = 0; i < N; i++)
                fprintf(fid, "%d %.9E\n", i, Mat[i]);
        fclose(fid);
}


bool to_bool(std::string str)
{
        std::transform(str.begin(), str.end(), str.begin(), ::tolower);
        std::istringstream is(str);
        bool b;
        is >> std::boolalpha >> b;
        return b;
}

void tokenparse(const string &input, string sep, vector<string> &tokens)
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