
#include "misc_func.hpp"

using namespace ::std;

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
        cout << "Reading int file " << filename << endl;
        FILE *fid;
        int n, m; //check size
        fid = fopen(filename.c_str(), "rb");
        fread(&n, sizeof(int), 1, fid);
        fread(&m, sizeof(int), 1, fid);
        if (n * m != Mat.size())
        {
                cout << "size does not match up resizing Mat to " << n * m << endl;
                Mat.resize(n * m, 0);
        }
        fread(Mat.data(), sizeof(int), n * m, fid);
        fclose(fid);
}

void writeMat(std::string filename, int m, int n, std::vector<int> &Mat)
{
        cout << "Writing int file " << filename << endl;
        FILE *fid;
        fid = fopen(filename.c_str(), "wb");
        fwrite(&n, sizeof(int), 1, fid);
        fwrite(&m, sizeof(int), 1, fid);
        if (n * m != Mat.size())
        {
                cout << "size does not match up in parameters " << n * m << endl;
                throw(-1);
        }
        fwrite(Mat.data(), sizeof(int), n * m, fid);
        fclose(fid);
}