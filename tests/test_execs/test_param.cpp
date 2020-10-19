

#include "param.hpp"
#include <assert.h>


using namespace::std;

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();
    if (rank != 0)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    paramMap p1("testInput.txt",rank);

    int testInt=p1.getParamInt("testInt");
    cout<<"test int is "<<testInt<<endl;
    assert(testInt==5);

    testInt=p1.getParamInt("testIntSpace");
    cout<<"test int Space is "<<testInt<<endl;
    assert(testInt==5);

    double testDouble=p1.getParamDouble("testDouble");
    cout<<"test double is "<<testDouble<<endl;
    assert(testDouble==3.59696);
    testInt=p1.getParamInt("testDoubleSpace");
    cout<<"test double Space is "<<testDouble<<endl;
    assert(testDouble==3.59696);

    string testString=p1.getParamString("testString");
    cout<<"test String is "<<testString<<endl;
    assert(testString==HelloWorld!);
    testString=p1.getParamString("testStringSpace");
    cout<<"test String Space is "<<testString<<endl;
    assert(testString==HelloWorld!);


    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}