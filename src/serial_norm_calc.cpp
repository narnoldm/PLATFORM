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
    string avgFile;
    if(argc>2)
    {
        debug_proc=atoi(argv[2]);
        if (argc > 3) {
            avgFile = argv[3];
        }

    }
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }
    string input = argv[1];
    cout<<"input string is: "<<input<<endl;
    vector<string> tokens;
    tokenparse(input,"|",tokens);
    tokens.push_back("Static_Pressure");
    tokens.push_back("-1");
    tokens.push_back("W");
    tokens.push_back("-1");
    tokens.push_back("Temperature");
    tokens.push_back("-1");
    tokens.push_back("Flamelet_Scalar_Mean");
    tokens.push_back("-1");
    tokens.push_back("Flamelet_Scalar_Variance");
    tokens.push_back("-1");
    tokens.push_back("Flamelet_Parameter");
    tokens.push_back("-1");


    tecIO *dataset1=new tecIO(tokens);
    //check for non standard input file name

    vector<double> data,average,norm,sum;
    int N=dataset1->nCells*dataset1->numVars;
    int M= dataset1->nSets;
    data.resize(N);
    average.resize(N);
    sum.resize(N);

    cout.precision(dbl::max_digits10);
    double val = 0.0;
    if (argc < 4) {
        for(int i=dataset1->snap0;i<=dataset1->snapF;i=i+dataset1->snapSkip)
        {
            dataset1->readSingle(i,data.data());
            
            for(int j=0;j<data.size();j++)
                average[j]+=data[j];
        }
        for(int j=0;j<data.size();j++)
            average[j]/=M;

        cout<<average[0]<<endl;

        for(int i=dataset1->snap0;i<=dataset1->snapF;i=i+dataset1->snapSkip)
        {
            dataset1->readSingle(i,data.data());
            for(int j=0;j<data.size();j++)
                sum[j]+=(data[j]-average[j])*(data[j]-average[j]);
        }
        cout<<sum[0]<<endl;
        for(int j=0;j<data.size();j++)
                sum[j]/=M;

        cout<<sum[0]<<endl;
        for(int j=0;j<data.size();j++)
            val+=sum[j];

	cout<<val<<endl;
        val/=N;
        val=sqrt(val);

    } else {
        // load average file
        dataset1->readAvg(avgFile);

        for(int i = dataset1->snap0; i <= dataset1->snapF; i= i + dataset1->snapSkip)
        {
            dataset1->readSingle(i,data.data());
            cout<<sum[0]<<endl;           
            cout<<data[0]<<"-"<<dataset1->average[0]<<"="<<data[0] - dataset1->average[0]<<endl;
            for(int j = 0; j < data.size(); j++) {
                sum[j] += (data[j] - dataset1->average[j])*(data[j] - dataset1->average[j]);
            }
        }
    	cout<<sum[0]<<endl;
        for(int j=0;j<data.size();j++)
                sum[j]/=M;
	cout<<sum[0]<<endl;
        
	for(int k =0; k<dataset1->numVars;k++)
	{
	val =0;
	for(int j = 0; j < dataset1->nCells; j++)
            val += sum[k*dataset1->nCells+j];

	cout<<val<<endl;
        val /= dataset1->nCells;
        val = sqrt(val);
	cout<<val<<endl;
	}
    }


    cout<<val<<endl;
    MPI_Barrier(MPI_COMM_WORLD);
    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}
