#include "metadata.hpp"
#include "param.hpp"

using namespace ::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    assert(argc == 3);

    int debug_proc = atoi(argv[1]);
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    string inpFile = argv[2];
    paramMap input(inpFile);

    string FOM;
    input.getParamString("FOM", FOM);
    string ROM;
    input.getParamString("ROM", ROM);

    string basis;
    input.getParamString("V", basis);

    string centerFile;
    input.getParamString("CenterFile", centerFile);
    string hashfile;
    input.getParamString("HashFile", hashfile);

    int FOMproj;
    input.getParamInt("FOM/ROM", FOMproj); //1 for pojected FOM 2 for ROM
    int comp;
    input.getParamInt("Comp", comp); //1 for normalized 2 for unnormalized 3 for uncentered

    int outPert;
    input.getParamInt("outPert", outPert);

    int outProj;
    input.getParamInt("outProj", outProj);
    int outErr;
    input.getParamInt("outErr", outErr);

    int subVar;
    input.getParamInt("SubVar", subVar); //0 for no 1 for yes
    //int outRecon=input.getParamInt("outRecon");

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);
    vector<string> token;

    cout << FOM << endl;
    tokenparse(FOM, "|", token);
    tecIO set1(token);

    if (FOMproj == 2)
        tokenparse(ROM, "|", token);

    tecIO set2(token);

    set1.activateReorder(hashfile);
    set2.activateReorder(hashfile);

    token.clear();

    tokenparse(basis, "|", token);
    meta SpaModes(token);
    token.clear();

    pMat V(SpaModes.nPoints, SpaModes.nSets, evenG, 0, 0, 0.0);

    SpaModes.batchRead(&V);
    //V.write_bin("V.bin");
    pMat q(SpaModes.nPoints, set1.nSets, evenG, 0, 0, 0.0);
    pMat VTq(SpaModes.nSets, set1.nSets, evenG, 0, 0, 0.0);
    pMat VVTq(SpaModes.nPoints, set1.nSets, evenG, 0, 0, 0.0);
    pMat pmVVTq(SpaModes.nPoints, set1.nSets, evenG, 0, 0, 0.0);

    //vector<double> err(set1.nSets,0.0);
    //vector<double> norm(set1.nSets,0.0);

    set1.calcCentering(&q, centerFile);
    string scaleMethod;
    input.getParamString("scaleMethod", scaleMethod);
    set1.calcScaling(&q, scaleMethod);

    set1.batchRead(&q);
    set1.centerData(&q);
    if (outPert == 1)
        set1.batchWrite(&q, "Pert", "Pert_");

    set1.scaleData(&q);

    cout << "computing VTq" << endl;
    VTq.matrix_Product('T', 'N', SpaModes.nSets, set1.nSets, SpaModes.nPoints, &V, 0, 0, &q, 0, 0, 1.0, 0.0, 0, 0);
    cout << "computing VVTq" << endl;
    VVTq.matrix_Product('N', 'N', SpaModes.nPoints, set1.nSets, SpaModes.nSets, &V, 0, 0, &VTq, 0, 0, 1.0, 0.0, 0, 0);
    cout << "done" << endl;

    if (comp >= 2)
    {
        set1.scaleData(&VVTq, true);
        if (comp >= 3)
            set1.centerData(&VVTq, true);
    }

    set2.batchRead(&q);
    if (comp < 3)
    {
        set1.centerData(&q);
        if (comp < 2)
        {
            set1.scaleData(&q);
        }
    }
    for (int k = 0; k < q.nelements; k++)
    {
        pmVVTq.dataD[k] = std::fabs(q.dataD[k] - VVTq.dataD[k]);
    }
    if (outProj)
        set1.batchWrite(&VVTq, "ProjectedSol_" + std::to_string(SpaModes.nSets), "Projeectedsol_");
    if (outErr)
        set1.batchWrite(&pmVVTq, "Error_" + std::to_string(SpaModes.nSets), "Error_");

    for (int k = 0; k < q.nelements; k++)
    {
        pmVVTq.dataD[k] = pmVVTq.dataD[k] * pmVVTq.dataD[k];
        q.dataD[k] = q.dataD[k] * q.dataD[k];
    }

    if (subVar == 0)
    {
        pMat ones(1, set1.nPoints, evenG, 0, 0, 1.0);
        pMat err(1, set1.nSets, evenG, 0, 0, 0.0);

        err.matrix_Product('N', 'N', 1, set1.nSets, set1.nPoints, &ones, 0, 0, &pmVVTq, 0, 0, 1.0, 0.0, 0, 0);
        for (int i = 0; i < err.nelements; i++)
        {
            err.dataD[i] = std::sqrt(err.dataD[i]);
        }
        err.write_bin("errorVector" + std::to_string(SpaModes.nSets) + ".bin");
    }
    if (subVar == 1)
    {
        pMat ones(1, set1.nCells, evenG, 0, 0, 1.0);
        pMat err(1, set1.nSets, evenG, 0, 0, 0.0);
        pMat norm(1, set1.nSets, evenG, 0, 0, 0.0);

        for (int v = 0; v < set1.numVars; v++)
        {
            err.matrix_Product('N', 'N', 1, set1.nSets, set1.nCells, &ones, 0, 0, &pmVVTq, v * set1.nCells, 0, 1.0, 0.0, 0, 0);
            if (FOMproj == 1)
                norm.matrix_Product('N', 'N', 1, set1.nSets, set1.nCells, &ones, 0, 0, &q, v * set1.nCells, 0, 1.0, 0.0, 0, 0);
            else if (FOMproj == 2)
                norm.matrix_Product('N', 'N', 1, set1.nSets, set1.nCells, &ones, 0, 0, &VVTq, v * set1.nCells, 0, 1.0, 0.0, 0, 0);

            for (int i = 0; i < err.nelements; i++)
            {
                norm.dataD[i] = std::sqrt(norm.dataD[i]);
                err.dataD[i] = std::sqrt(err.dataD[i]) / norm.dataD[i];
            }
            err.write_bin("error_" + set1.varName[v] + std::to_string(SpaModes.nSets) + ".bin");
        }
    }

    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}