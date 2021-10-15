
#include "misc_func.hpp"
#include <cmath>

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

/**
 * Computes aggregate timing statistics across all MPI processes. 
 * 
 * Calculates average, minimum, maximum, and standard deviation.
 * Appends recorded results to end of log file given by filename, with label given by label
*/
void aggregateTiming(double t, string filename, string label) {

	int rank, num_procs;
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// average
	double average;
	MPI_Allreduce(&t, &average, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD);
	average /= num_procs;
     
    // min and max
	double maximum, minimum;
    MPI_Reduce(&t, &maximum, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t, &minimum, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD);

    // standard deviation
	double global_squareDiff;
    double local_squareDiff = pow(t - average, 2.0);
    MPI_Reduce(&local_squareDiff, &global_squareDiff, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double standard_dev = sqrt(global_squareDiff / double(num_procs - 1));
		ofstream out;
		out.open(filename, ios::app);
		out << label << " timings" << endl;
		out << "==================" << endl;
		out << "Average: " << to_string(average) << endl;
		out << "Maximum: " << to_string(maximum) << endl;
		out << "Minimum: " << to_string(minimum) << endl;
		out << "STD: " << to_string(standard_dev) << endl;
		out.close();
	}

}