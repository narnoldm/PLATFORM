

#include "parser.hpp"

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::ofstream sink("/dev/null");
    streambuf* strm_buffer = cout.rdbuf();
    if(rank!=0)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    cout << "initializing" << endl;

    inputReader* file1p= new inputReader("input.pdp");
    executioner *exec = new executioner(file1p);
    
    exec->init();
    exec->exec_all();
    exec->output();
    exec->clear();



    
    delete file1p;
    MPI_Barrier(MPI_COMM_WORLD);
    cout.rdbuf (strm_buffer);
    MPI_Finalize();
    return 0;
}