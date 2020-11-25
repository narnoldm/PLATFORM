

#include "pMat.hpp"
#include "metadata.hpp"
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

    if(argc<3)
    {
        MPI_Abort(MPI_COMM_WORLD,-1);
    }
    int debug_proc=atoi(argv[1]);
    int M=atoi(argv[2]);
    int N=atoi(argv[3]);


    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    cout<<" M "<<M<<endl<<" N "<<N<<endl;

    if((M*N)>1e6)
    {
        cout<< (M*N)/(1e6)*8<< " MB"<<endl;
        cout<<"Asking for very large matrix are you sure you are on a compute node"<<endl;
        if(rank==debug_proc)
        {
            cout<<"hit enter to continue"<<endl;
            cin.get();
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    PGrid *p1 = new PGrid(rank,size,0);

    pMat *A=new pMat(M,N,p1);

    for(int i=0;i<A->nelements;i++)
    {
        A->dataD[i]=rand();
    }
    A->write_bin("A.bin");
    A->read_bin("A.bin");

    pMat *U,*VT;
    U=new pMat(A->M,min(A->M,A->N),p1); 
    VT=new pMat(min(A->M,A->N),A->N,p1); 
    vector<double> S(min(A->M,A->N));

    meta *m1,*m2;
    string prefix="a";
    string suffix=".bin";
    m1= new meta();
    m1->snap0 = 0;
    m1->snapF = A->N;
    m1->snapSkip = 1;
    m1->prefix = prefix;
    m1->suffix = suffix;
    m1->nPoints = A->M;
    m1->nSets=A->N;
    m1->batchWrite(A);

    string prefix2="out/a";
    string suffix2=".bin";
    m2=new meta(0,A->N-1,1,prefix2,suffix2);
    m2->batchRead(A);



    A->svd_run(A->M,A->N,0,0,U,VT,S);




    delete A;
    delete m1;
    delete m2;
    delete p1;

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}