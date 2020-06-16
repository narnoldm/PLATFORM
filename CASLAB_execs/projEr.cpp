#include "metadata.hpp"
#include "param.hpp"


using namespace::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    assert(argc==3);

    int debug_proc=atoi(argv[1]);
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }
    
    string inpFile=argv[2];
    paramMap input(inpFile);

    string FOM = input.getParamString("FOM");
    string ROM = input.getParamString("ROM");

    string basis=input.getParamString("V");

    string centering=input.getParamString("CenterFile");
    string hashfile=input.getParamString("HashFile");

    PGrid *evenG;
    evenG = new PGrid(rank,size,0);
    vector<string> token;




    
    tokenparse(FOM,"|",token);
    tecIO set1(token);
    
    tokenparse(ROM,"|",token);
    tecIO set2(token);
    
    set1.activateReorder(hashfile);
    set2.activateReorder(hashfile);

    token.clear();


    tokenparse(basis,"|",token);
    meta SpaModes(token);
    token.clear();

    pMat V(SpaModes.nPoints,SpaModes.nSets,evenG,0,0,0.0);

    SpaModes.batchRead(&V); 


    pMat q(SpaModes.nPoints,1,evenG,0,0,0.0);
    pMat VTq(SpaModes.nSets,1,evenG,0,0,0.0);
    pMat VVTq(SpaModes.nPoints,1,evenG,0,0,0.0);
    

    vector<double> err(set1.nSets,0.0);

    set1.readAvg(centering);
    set1.calcNorm(&q);

    for(int i=0;i<set1.nSets;i++)
    {
        set1.batchRead(&q,i);
        set1.normalize(&q); 
        VTq.matrix_vec_product('T',SpaModes.nSets,SpaModes.nPoints,1.0,&V,0,0,&q,0,0,0.0,0,0);
        VVTq.matrix_vec_product('N',SpaModes.nPoints,SpaModes.nSets,1.0,&V,0,0,&VTq,0,0,0.0,0,0);
        set2.batchRead(&q,i);
        set1.normalize(&q);
        for(int k=0;k<q.nelements;k++)
        {
            q.dataD[k]-=VVTq.dataD[k];
            err[i]+=std::pow(q.dataD[k],2); 
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,err.data(),err.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(int i=0;i<err.size();i++)
        err[i]=std::sqrt(err[i]);

    if(rank==0)
        printASCIIVecP0("projErr.txt",err.data(),err.size());

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0; 
}