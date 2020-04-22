#include "inputReader.hpp"

#include <limits>
typedef numeric_limits<double> dbl;

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
    string input = argv[1];
    cout<<"input string is: "<<input<<endl;
    vector<string> tokens;
    tokenparse(input,"|",tokens);

    tecIO *dataset1=new tecIO(tokens);
    //check for non standard input file name

    vector<double> data,average,norm,sum;
    int N=dataset1->nCells*dataset1->numVars;
    int M= dataset1->nSets;
    data.resize(N);
    average.resize(N);
    sum.resize(N);

    for(int i=dataset1->snap0;i<=dataset1->snapF;i=i+dataset1->snapSkip)
    {
        dataset1->readSingle(i,data.data());
        for(int j=0;j<data.size();j++)
            average[j]+=data[j];
    }
    for(int j=0;j<data.size();j++)
            average[j]/=M;

    cout.precision(dbl::max_digits10);
    cout<<average[0]<<endl;

    for(int i=dataset1->snap0;i<=dataset1->snapF;i=i+dataset1->snapSkip)
    {
        dataset1->readSingle(i,data.data());
        for(int j=0;j<data.size();j++)
            sum[j]+=(data[j]-average[j])*(data[j]-average[j]);
    }
    for(int j=0;j<data.size();j++)
            sum[j]/=M;

    double val=0.0;
    for(int j=0;j<data.size();j++)
        val+=sum[j];

    val/=N;
    val=sqrt(val);


    cout<<val<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}