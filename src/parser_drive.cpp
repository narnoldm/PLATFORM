

#include "parser.hpp"

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    cout << "initializing" << endl;

    inputReader file1("input.pdp");
    inputReader* file1p= &file1;
 
    executioner exec(file1p);
    exec.init();
    exec.exec_all();
    exec.clear();

    MPI_Finalize();
    return 0;
}