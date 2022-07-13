#include "processGrid.hpp"

using namespace :: std;

PGrid::PGrid(int r, int s, int type)
{
    rank = r;
    size = s;
    printRank = (rank == 0);
    pdims[0] = 0;
    pdims[1] = 0;
    icntxt = 0;
    myrow = 0;
    mycol = 0;
    int ngone = -1, zero = 0;
    string orderStr = "R";
    char *order = new char[orderStr.size() + 1];
    copy(orderStr.begin(), orderStr.end(), order);
    order[orderStr.size()] = '\0';

    // dimensions of process grid
    if (type == 0)
    {
        MPI_Dims_create(size, 2, pdims);
    }
    else if (type == 1)
    {
        pdims[0] = size;
        pdims[1] = 1;
    }
    else if (type == 2)
    {
        pdims[0] = 1;
        pdims[1] = 1;
    }
    else if (type == 3)
    {
        pdims[0] = 1;
        pdims[1] = size;
    }

    cout << "Initializing Cblacs" << endl;
    Cblacs_pinfo(&rank, &size);
    Cblacs_get(-1, 0, &icntxt);
    Cblacs_gridinit(&icntxt, order, pdims[0], pdims[1]);
    Cblacs_gridinfo(icntxt, &pcol, &prow, &myrow, &mycol);

    cout << "Processor Grid M(cols)=" << pdims[0] << " N(rows)=" << pdims[1] << endl;
    cout << "local processor row " << myrow << " and column " << mycol << endl;
    cout << "PBLACS is row major by default and everything is else is Column major (because reasons)" << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    delete[] order;
}
PGrid::PGrid(int r, int s)
{



}
PGrid::~PGrid()
{
        cout << "clearing pblacs buffers" << endl;
        Cblacs_gridexit(icntxt);
}
int PGrid::getDim(int dim)
{
        return pdims[dim];
}