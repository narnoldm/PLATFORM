
#include "misc_func.hpp"

using namespace ::std;

void printASCIIVecP0(std::string fname, double *Mat, int N)
{
        FILE *fid;
        fid = fopen(fname.c_str(), "w");
        for (int i = 0; i < N; i++)
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

void readMat(std::string filename, std::vector<int> &Mat)
{
        printf("Reading int file %s\n", filename.c_str());
        FILE *fid;
        int n, m; //check size
        fid = fopen(filename.c_str(), "rb");
        fread(&n, sizeof(int), 1, fid);
        fread(&m, sizeof(int), 1, fid);
        if (n*m != Mat.size())
        {
                printf("size does not match up resizing Mat\n");
                Mat.resize(n*m,0);
        }
        fread(Mat.data(),sizeof(int),n*m,fid);
        fclose(fid);
}