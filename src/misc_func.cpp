
#include "misc_func.hpp"
#include <cmath>

using namespace ::std;

// write array of doubles (of at least length N) to file
void printASCIIVecP0(std::string fname, double *Mat, int N)
{
    FILE *fid;
    fid = fopen(fname.c_str(), "w");
    for (int i = 0; i < N; i++)
        fprintf(fid, "%d %.9E\n", i, Mat[i]);
    fclose(fid);
}


// convert string to boolean
// handles any case differences in string, maps to "true" or "false"
bool to_bool(std::string str)
{
    std::transform(str.begin(), str.end(), str.begin(), ::tolower);
    std::istringstream is(str);
    bool b;
    is >> std::boolalpha >> b;
    return b;
}

// break string into a vector of separate "tokens", delimited by sep
void tokenparse(const string &input, string sep, vector<string> &tokens)
{
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

// read vector of ints from binary file
// File should have two 4-byte int header values defining the number of rows and columns
// of the matrix which is defined by the vector
void readMat(std::string filename, std::vector<int> &Mat)
{
    cout << "Reading int file " << filename << endl;
    FILE *fid;
    int n, m; //check size
    fid = fopen(filename.c_str(), "rb");
    size_t warn;
    warn = fread(&n, sizeof(int), 1, fid);
    warn = fread(&m, sizeof(int), 1, fid);
    if (n * m != Mat.size())
    {
        cout << "size does not match up resizing Mat to " << n * m << endl;
        Mat.resize(n * m, 0);
    }
    warn = fread(Mat.data(), sizeof(int), n * m, fid);
    fclose(fid);
}

// write vector of ints to binary file,
// with two 4-byte int header values defining the number of rows and columns
// of the matrix which is defined by the vector
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

// write vector of doubles to ASCII file, with single string header
void writeASCIIDoubleVec(std::string filename, std::vector<double> &vec)
{
    FILE *fid;
    if ((fid = fopen(filename.c_str(), "w")) == NULL)
    {
        printf("error with file open\n");
    }
    fprintf(fid, "file=%s\n", filename.c_str());
    for (int i = 0; i < vec.size(); ++i)
    {
        fprintf(fid, "%16.16E\n", vec.data()[i]);
    }
    fclose(fid);
}
