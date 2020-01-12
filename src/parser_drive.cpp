#include <iostream>

//#include "parser.cpp"
#include "parser.hpp"
using namespace ::std;

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "initializing" << endl;

    string teststring = "input|tecplot|gemsma_cmb_|.szplt|1|100|1";
    string sep = "|";

    vector<string> tokenlist;

    tokenparse(teststring, sep, tokenlist);

    for (int i = 0; i < tokenlist.size(); i++)
        cout << tokenlist[i] << endl;

    inputReader file1("input.pdp");

    inputReader* file1p= &file1;

    executioner exec(file1p);
    exec.init();
    exec.exec_all();
    exec.clear();

    MPI_Finalize();
    return 0;
}