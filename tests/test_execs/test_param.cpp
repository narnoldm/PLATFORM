

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

    int testInt;
    assert(p1.getParamInt("testInt",testInt));
    cout<<"test int is "<<testInt<<endl;
    assert(testInt==5);

    assert(p1.getParamInt("testIntSpace",testInt));
    cout<<"test int Space is "<<testInt<<endl;
    assert(testInt==5);

    double testDouble;
    assert(p1.getParamDouble("testDouble",testDouble));
    cout<<"test double is "<<testDouble<<endl;
    assert(testDouble==3.59696);
    assert(p1.getParamDouble("testDoubleSpace",testDouble));
    cout<<"test double Space is "<<testDouble<<endl;
    assert(testDouble==3.59696);

    string testString;
    assert(p1.getParamString("testString",testString));
    cout<<"test String is "<<testString<<endl;
    assert(testString=="HelloWorld!");
    assert(p1.getParamString("testStringSpace",testString));
    cout<<"test String Space is "<<testString<<endl;
    assert(testString=="HelloWorld!");

    bool testBool;
    assert(p1.getParamBool("testBool",testBool));
    cout<<"test Bool is "<<testBool<<endl;
    assert(testBool==true);
    assert(p1.getParamBool("testBoolSpace",testBool));
    cout<<"test Bool Space is "<<testBool<<endl;
    assert(testBool==true);

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}