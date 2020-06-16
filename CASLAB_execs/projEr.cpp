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

    //int outRecon=input.getParamInt("outRecon");

    PGrid *evenG;
    evenG = new PGrid(rank,size,0);
    vector<string> token;




    cout<<FOM<<endl;
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
    //V.write_bin("V.bin");
    pMat q(SpaModes.nPoints,1,evenG,0,0,0.0);
    pMat VTq(SpaModes.nSets,1,evenG,0,0,0.0);
    pMat VVTq(SpaModes.nPoints,1,evenG,0,0,0.0);
    

    vector<double> err(set1.nSets,0.0);

    set1.readAvg(centering);
    set1.calcNorm(&q);

    //set1.activateGEMSbin(hashfile);
    //set1.meshFile=hashfile;
    //set1.fixedMesh=true;
    for(int i=0;i<set1.nSets;i++)
    {
        set1.batchRead(&q,i);
        set1.subAvg(&q);
        set1.normalize(&q); 
        //q.write_bin("q"+std::to_string(i)+".bin");
        VTq.matrix_Product('T','N',SpaModes.nSets,1,SpaModes.nPoints,&V,0,0,&q,0,0,1.0,0.0,0,0);
        //VTq.write_bin("VTq"+std::to_string(i)+".bin");
        //V.write_bin("Vafm"+std::to_string(i)+".bin");
        VVTq.matrix_Product('N','N',SpaModes.nPoints,1,SpaModes.nSets,&V,0,0,&VTq,0,0,1.0,0.0,0,0);
        //VVTq.write_bin("VVTq"+std::to_string(i)+".bin");
        set2.batchRead(&q,i);
        set1.subAvg(&q);
        set1.normalize(&q);
        for(int k=0;k<q.nelements;k++)
        {
            q.dataD[k]-=VVTq.dataD[k];
            /*cout<<VTq.dataD[k]<<endl;
            cout<<V.dataD[k]<<endl;
            cout<<VVTq.dataD[k]<<endl;*/
            err[i]+=q.dataD[k]*q.dataD[k]; 
        }
    }
    MPI_Allreduce(MPI_IN_PLACE,err.data(),err.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(int i=0;i<err.size();i++)
        err[i]=std::sqrt(err[i]);//)/set1.nPoints;

    if(rank==0)
        printASCIIVecP0("projErr"+std::to_string(SpaModes.nSets)+".txt",err.data(),err.size());
    for(int i=0;i<err.size();i++)
        err[i]=err[i]/set1.nPoints;
    printASCIIVecP0("projErrN"+std::to_string(SpaModes.nSets)+".txt",err.data(),err.size());
    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0; 
}