

#include "parser.hpp"

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    //silence output for non root MPI processes
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();
    int debug_proc=0;
    if(argc>2)
    {
        debug_proc=atoi(argv[2]);
    }
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    //check for non standard input file name
    string input;
    cout<<argc;
    input="input.pdp";
    if(argc>1)
    {
        string input = argv[1];
    }
    cout << input << endl;
    cout << "initializing" << endl;



    //launch parser
    inputReader *file1p = new inputReader(input);
    executioner *exec = new executioner(file1p);

    exec->init();
    //exec->exec_all();
    //exec->output();
    exec->clear();



    delete file1p,exec;
    MPI_Barrier(MPI_COMM_WORLD);
    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}