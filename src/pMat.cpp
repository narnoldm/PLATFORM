#ifndef P_MAT_H
#define P_MAT_H
#include "pMat.hpp"
#endif

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
        char *order = "R";
        if (type == 0)
                MPI_Dims_create(size, 2, pdims);
        else if (type == 1)
        {
                pdims[0] = 1;
                pdims[1] = size;
        }
        else if (type == 2)
        {
                pdims[0] = 1;
                pdims[1] = 1;
        }
        if (printRank)
                printf("Processor Grid N=%d M=%d\n", pdims[0], pdims[1]);
        Cblacs_pinfo(&rank, &size);
        Cblacs_get(-1, 0, &icntxt);
        Cblacs_gridinit(&icntxt, order, pdims[0], pdims[1]);
        Cblacs_gridinfo(icntxt, &prow, &pcol, &myrow, &mycol);
        MPI_Barrier(MPI_COMM_WORLD);
}
PGrid::~PGrid()
{
        if (printRank)
                printf("clearing pblacs buffers\n");
        Cblacs_gridexit(icntxt);
}
int PGrid::getDim(int dim)
{
        return pdims[dim];
}
pMat::pMat()
{
        cout << "empty pMat created" << endl;
}
pMat::pMat(int n, int m, PGrid *pGp, int t, int b, double init)
{
        pG = pGp;
        setupMat(n, m, t, b, 1, init);
}
pMat::pMat(int n, int m, PGrid *pGp, int t, int b, int c, double init)
{
        pG = pGp;
        setupMat(n, m, t, b, c, init);
}
pMat::~pMat()
{
        if (printRank)
                printf("Deallocating Distributed Matrix\n");
}
void pMat::setupMat(int n, int m, int t, int b, int c, double init)
{
        N = n;
        M = m;
        type = t;
        block = b;
        cycles = c;
        printRank = pG->printRank;

        if (printRank)
                printf("Creating Matrix\nN=%d M=%d\n", N, M);
        if (block == 0) //square blocks
        {
                nb = 1000000;
                mb = nb;
                if (nb > (N / pG->getDim(0)))
                        nb = std::max(1, N / (pG->getDim(0) * cycles));
                if (nb > (M / pG->getDim(1)))
                        nb = std::max(1, M / (pG->getDim(1) * cycles));
                mb = nb;
                if (pG->printRank)
                {
                        printf("mb/nb = %d\n", nb);
                }
                MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (block == 1) //other blocks
        {
                nb = 1000000;
                mb = nb;
                if (nb > (N / pG->getDim(0)))
                        nb = N / pG->getDim(0);
                if (mb > (M / pG->getDim(1)))
                        mb = M / pG->getDim(1);
                if (pG->printRank)
                {
                        printf("nb is %d\n", nb);
                        printf("mb is %d\n", mb);
                }
                MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (block == 2) //p0 block
        {
                nb = N;
                mb = M;
                if (pG->printRank)
                {
                        printf("nb is %d\n", nb);
                        printf("mb is %d\n", mb);
                        fflush(stdout);
                }
                MPI_Barrier(MPI_COMM_WORLD);
        }
        myRC[0] = numroc_(&N, &nb, &(pG->myrow), &i_zero, &(pG->prow));
        myRC[1] = numroc_(&M, &mb, &(pG->mycol), &i_zero, &(pG->pcol));
        if (myRC[0] < 1)
                myRC[0] = 1;
        if (myRC[1] < 1)
                myRC[1] = 1;
        nelements = (long long)myRC[0] * (long long)myRC[1];
        if (type == 0)
        {
                if (printRank)
                        printf("Mat is Double\n");
                dataD.resize(nelements, init);
                GBs=nelemets*8/(1e9);
        }
        else if (type == 1)
        {
                if (printRank)
                        printf("Mat is Complex\n");
                dataC.resize(nelements);
                for (int i = 0; i < nelements; i++)
                {
                        dataC[i].real = init;
                        dataC[i].imag = 0.0;
                        GBs=nelemets*16/(1e9);
                }
                GBs=nelments*
        }
        else
        {
                if(printRank)
                {
                        printf("invalid type\n");
                        throw(-1);
                }
        }
        
        
        descinit_(desc, &N, &M, &nb, &mb, &i_zero, &i_zero, &(pG->icntxt), &myRC[0], &info);
        if (info != 0)
        {
                if (pG->printRank)
                        printf("Error in descriptor setup in arguement %d\n", -info);
        }
        if (printRank)
                printf("Matrix Constructed\n");
}

void pMat::switchType(int t)
{
        if (type == t)
        {
                if (pG->printRank)
                        printf("matrix is already of data type %d\n", type);
                return;
        }
        else if (t == 0)
        {
                if (pG->printRank)
                        printf("Switching to Double\n");
                dataD.resize(nelements);
                for (int i = 0; i < nelements; i++)
                        dataD[i] = dataC[i].real;
                dataC.clear();
                type = 0;
                if (pG->printRank)
                        printf("Matrix is now type %d\n", type);
                return;
        }
        else if (t == 1)
        {
                if (pG->printRank)
                        printf("Switching to Complex\n");

                dataC.resize(nelements);
                for (int i = 0; i < nelements; i++)
                {
                        dataC[i].real = dataD[i];
                        dataC[i].imag = 0;
                }
                dataD.clear();
                type = 1;
                if (pG->printRank)
                        printf("Matrix is now type %d\n", type);
                return;
        }
}

void pMat::printMat()
{
        for (int p = 0; p < 1; p++)
        {
                if (p == pG->rank)
                {
                        printf("processor %d has :", p);
                        for (int i = 0; i < nelements; i++)
                        {
                                if (i % myRC[0] == 0)
                                        printf("\n");

                                if (type == 0)
                                        printf("%f\n", dataD[i]);
                                else if (type == 1)
                                        printf("%f + %f i\n", dataC[i].real, dataC[i].imag);
                        }
                        fflush(stdout);
                }
                MPI_Barrier(MPI_COMM_WORLD);
        }
}

int pMat::write_bin(std::string filename)
{
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File fH;
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fH);
        if (printRank)
        {
                MPI_File_write(fH, &N, 1, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_write(fH, &M, 1, MPI_INT, MPI_STATUS_IGNORE);
                printf("Write Start\n");
                printf("N=%d nb=%d, M=%d mb=%d\n", N, nb, M, mb);
        }
        MPI_File_close(&fH);
        MPI_Barrier(MPI_COMM_WORLD);
        int disp = 2 * sizeof(int);
        int dims[2] = {N, M};
        int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
        int dargs[2] = {nb, mb};
        MPI_Datatype darray;
        MPI_Type_create_darray(pG->size, pG->rank, 2, dims, distribs, dargs, pG->pdims, MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
        MPI_Type_commit(&darray);
        int tsize, mpiEls;
        MPI_Type_size(darray, &tsize);
        mpiEls = tsize / (sizeof(double));
        if (nelements != mpiEls)
        {
                printf("Allocation via MPI %d and pblacs %d is different\n", mpiEls, nelements);
                throw(-1);
        }
        double t2, t1;
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fH);
        if (printRank)
                printf("MPI Allocation %d , pblacs Allocation %d\n", mpiEls, nelements);
        t1 = MPI_Wtime();
        MPI_File_set_view(fH, disp, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
        MPI_File_write_all(fH, dataD.data(), mpiEls, MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_close(&fH);
        t2 = MPI_Wtime();
        if (printRank)
                printf("Write time is %f\n", t2 - t1);
}


int pMat::read_bins(std::string prefix, int start, int skip, int end)
{
        //checks
        if (pG->prow != 1)
        {
                printf("blocking must be rectangular for benign read in\n");
                throw(-1);
        }
        int numFiles = nelements / N;
        int iP = 0;
        int fileIndex = 0;
        int localC = 0;
        std::ostringstream name;
        for (int i = 0; i < M; i++)
        {
                iP = (int)(i / mb);
                while (iP > (pG->size - 1))
                {
                        iP = iP - pG->size;
                }
                if (pG->rank == iP)
                {
                        name.str("");
                        name.clear();
                        name << std::fixed;
                        name << std::setprecision(6);
                        fileIndex = start + i * skip;
                        if (printRank)
                                printf("proc %d is reading file %d\n", iP, fileIndex);
                        name << fileIndex;
                        read_single_bin((prefix + name.str() + ".bin").c_str(), localC);
                        localC++;
                }
        }
}

int pMat::read_single_bin(const char *name, int col)
{
        FILE *fid;
        fid = fopen(name, "rb");
        int Mfile, Nfile;
        fread(&Mfile, sizeof(int), 1, fid);
        fread(&Nfile, sizeof(int), 1, fid);
        if (Mfile != 0)
        {
                printf("file has more than 1 snapshot !!!\n");
                return -1;
        }
        if (Nfile != N)
        {
                printf("file has incorrect number of points");
                return -1;
        }
        double *pointD = dataD.data() + N * col;
        fread(pointD, sizeof(double), N, fid);
        fclose(fid);
        return 1;
}

int pMat::read_bin(string &filename)
{

        int rN, rM;
        double t2, t1;
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File fH;
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fH);
        MPI_File_read_all(fH, &rN, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read_all(fH, &rM, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fH);

        while ((rN != N) || (rM != M))
        {
                if (printRank)
                        printf("bin file and matrix do not have same dimension\n Attempting to reassign Matrix\n");

                this->~pMat();
                this->setupMat(rN, rM, type, block, cycles, 0.0);
        }
        int dims[2] = {N, M};
        int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
        int dargs[2] = {nb, mb};
        MPI_Datatype darray;
        MPI_Type_create_darray(pG->size, pG->rank, 2, dims, distribs, dargs, pG->pdims, MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
        MPI_Type_commit(&darray);
        int disp = 2 * sizeof(int);
        int tsize, mpiEls;
        MPI_Type_size(darray, &tsize);
        mpiEls = tsize / (sizeof(double));
        //check size
        if (nelements != mpiEls)
        {
                printf("Allocation via MPI %d and pblacs %d is different\n", mpiEls, nelements);
                throw(-1);
        }
        if (printRank)
                printf("MPI Allocation %d , pblacs Allocation %d\n", mpiEls, nelements);
        if (printRank)
                printf("Read Starting\n");
        t1 = MPI_Wtime();
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fH);
        MPI_File_set_view(fH, disp, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
        MPI_File_read_all(fH, dataD.data(), mpiEls, MPI_DOUBLE, MPI_STATUS_IGNORE);
        MPI_File_close(&fH);
        t2 = MPI_Wtime();
        if (printRank)
                printf("read of %s took %f seconds\n", filename.c_str(), t2 - t1);

        return 1;
}

bool pMat::check_bin_size(string filename, int &mN, int &mM)
{
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File fH;
        if (!MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fH))
        {
                MPI_File_read_all(fH, &mN, 1, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_read_all(fH, &mM, 1, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_close(&fH);
                cout << filename << " has data with N: " << mN << " M: " << mM << endl;
                return true;
        }
        else
        {
                cout << filename << " not found" << endl;
                return false;
        }
}

int pMat::matrix_Product(char tA, char tB, int n, int m, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc)
{

        if ((A->type == 0) && (B->type == 0) && (type == 0))
        {
                //if(printRank)
                //        printf("Double Multiply\n");
                int IA = ia + 1;
                int JA = ja + 1;
                int IB = ib + 1;
                int JB = jb + 1;
                int IC = ic + 1;
                int JC = jc + 1;
                pdgemm(&tA, &tB, &n, &m, &k, &alpha, A->dataD.data(), &IA, &JA, A->desc, B->dataD.data(), &IB, &JB, B->desc, &beta, dataD.data(), &IC, &JC, desc);
        }
        else if ((A->type == 1) && (B->type == 1) && (type == 1))
        {
                if (printRank)
                        printf("Complex Multiply\n");
                int IA = ia + 1;
                int JA = ja + 1;
                int IB = ib + 1;
                int JB = jb + 1;
                int IC = ic + 1;
                int JC = jc + 1;
                MKL_Complex16 Ac, Bc;
                Ac.real = alpha;
                Ac.imag = 0;
                Bc.real = beta;
                Bc.imag = 0;
                pzgemm(&tA, &tB, &n, &m, &k, &Ac, A->dataC.data(), &IA, &JA, A->desc, B->dataC.data(), &IB, &JB, B->desc, &Bc, dataC.data(), &IC, &JC, desc);
        }
        else
        {
                if (printRank)
                        printf("Other Formats not supported yet\n");
        }
        return 0;
}

int pMat::matrix_Sum(char tA, int n, int m, pMat *A, int ia, int ja, int ib, int jb, double alpha, double beta)
{

        if ((A->type == 0) && (type == 0))
        {
                if (printRank)
                        printf("Double Sum\n");
                int IA = ia + 1;
                int JA = ja + 1;
                int IB = ib + 1;
                int JB = jb + 1;
                pdgeadd(&tA, &n, &m, &alpha, A->dataD.data(), &IA, &JA, A->desc, &beta, dataD.data(), &IB, &JB, desc);
        }
        else if ((A->type == 1) && (type == 1))
        {
                if (printRank)
                        printf("Complex Sum\n");
                int IA = ia + 1;
                int JA = ja + 1;
                int IB = ib + 1;
                int JB = jb + 1;
                MKL_Complex16 Ac, Bc;
                Ac.real = alpha;
                Ac.imag = 0;
                Bc.real = beta;
                Bc.imag = 0;
                pzgeadd(&tA, &n, &m, &Ac, A->dataC.data(), &IA, &JA, A->desc, &Bc, dataC.data(), &IB, &JB, desc);
        }
        else
        {
                if (printRank)
                        printf("Other Formats not supported yet\n");
        }
        return 0;
}

int pMat::svd_run(int N, int M, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S)
{
        int info;
        char *JOBU = "V", *JOBVT = "V";
        std::vector<double> WORK(1);
        //double *WORK=(double*)mkl_calloc(1,sizeof(double),64);
        int LWORK = -1;
        int IA = ia + 1;
        int JA = ja + 1;
        int i_one = 1;
        double t2, t1;
        pdgesvd(JOBU, JOBVT, &N, &M, dataD.data(), &IA, &JA, desc, S.data(), U->dataD.data(), &IA, &JA, U->desc, VT->dataD.data(), &i_one, &i_one, VT->desc, WORK.data(), &LWORK, &info);
        if (printRank)
                printf("WORK=%f, LWORK=%d,info=%d\n", WORK[0], LWORK, info);
        LWORK = WORK[0];
        //mkl_free(WORK);
        WORK.resize(LWORK);
        //WORK=(double*)mkl_calloc(LWORK,sizeof(double),64);
        if (printRank)
                printf("Work Allocated: %f MB per processor\n", LWORK / (1e6) * 8);
        //SVD run
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        pdgesvd(JOBU, JOBVT, &N, &M, dataD.data(), &IA, &JA, desc, S.data(), U->dataD.data(), &IA, &JA, U->desc, VT->dataD.data(), &i_one, &i_one, VT->desc, WORK.data(), &LWORK, &info);
        t2 = MPI_Wtime();
        if (printRank)
                printf("SVD complete in %f seconds\n", t2 - t1);
        //mkl_free(WORK);
        return 1;
}
int pMat::transpose(pMat *A, int n, int m, int ia, int ja)
{
        int IA = ia + 1;
        int JA = ja + 1;
        if (printRank)
                printf("Copying transpose\n");

        if ((n != N) && (m != A->N))
                if (printRank)
                {
                        printf("transpose dimension mismatch\n");
                        return -1;
                }
        int i_one = 1;
        double ONE = 1.0;
        double ZERO = 0.0;
        pdtran(&n, &m, &ONE, A->dataD.data(), &IA, &JA, A->desc, &ZERO, dataD.data(), &i_one, &i_one, desc);
}
int pMat::changeContext(pMat *A, int n, int m, int ia, int ja, int ib, int jb)
{
        int IA = ia + 1;
        int JA = ja + 1;
        int IB = ib + 1;
        int JB = jb + 1;
        if (printRank)
                printf("Copying Matrix\n n=%d , m =%d\n", n, m);
        int i_one = 1;
        if (type == 0)
        {
                //if(printRank)
                //       printf("Double\n");
                pdgemr2d(&n, &m, A->dataD.data(), &IA, &JA, A->desc, dataD.data(), &IB, &JB, desc, &(pG->icntxt));
        }
        if (type == 1)
        {
                if (printRank)
                        printf("Complex\n");
                pzgemr2d(&n, &m, A->dataC.data(), &IA, &JA, A->desc, dataC.data(), &IB, &JB, desc, &(pG->icntxt));
        }
}

int pMat::changeContext(pMat *A)
{
        changeContext(A, N, M, 0, 0, 0, 0);
}
int pMat::dMax(int dim, int rc, double val)
{
        double max;
        int index = 0;

        if (printRank)
                printf("finding max\n");
        if (dim == 0)
        {
                int IA = 1, JA = rc + 1, i_one = 1;
                pdamax(&N, &val, &index, dataD.data(), &IA, &JA, desc, &i_one);
        }
        if (dim == 1)
        {
                int IA = rc + 1, JA = 1, i_one = 1;
                pdamax(&M, &val, &index, dataD.data(), &IA, &JA, desc, &N);
        }

        if (printRank)
                printf("max is %f at %d\n", val, index);
}
int pMat::dAve(int dim, int rc, double val)
{

        if (printRank)
                printf("finding sum\n");
        if (dim == 0)
        {
                int IA = 1, JA = rc + 1, i_one = 1;
                pdasum(&N, &val, dataD.data(), &IA, &JA, desc, &i_one);
        }
        if (dim == 1)
        {
                int IA = rc + 1, JA = 1, i_one = 1;
                pdasum(&M, &val, dataD.data(), &IA, &JA, desc, &N);
        }

        if (printRank)
                printf("sum is %f\n", val);
}

void pMat::printDesc()
{
        if (printRank)
        {
                std::cout << "Descriptor type: " << desc[0] << std::endl;
                std::cout << "BLACS context: " << desc[1] << std::endl;
                std::cout << "Global Rows: " << desc[2] << std::endl;
                std::cout << "Global Cols: " << desc[3] << std::endl;
                std::cout << "Row Blocking factor: " << desc[4] << std::endl;
                std::cout << "Column Blocking factor: " << desc[5] << std::endl;
                std::cout << "Process row where first row is: " << desc[6] << std::endl;
                std::cout << "Process Col where first col is: " << desc[7] << std::endl;
                std::cout << "Leading Dimension: " << desc[8] << std::endl;
        }
}
ostream &operator<<(std::ostream &os, const pMat &p)
{
        if(p.printRank)
        {
                std::cout << "Descriptor type: " << desc[0] << std::endl;
                std::cout << "BLACS context: " << desc[1] << std::endl;
                std::cout << "Global Rows: " << desc[2] << std::endl;
                std::cout << "Global Cols: " << desc[3] << std::endl;
                std::cout << "Row Blocking factor: " << desc[4] << std::endl;
                std::cout << "Column Blocking factor: " << desc[5] << std::endl;
                std::cout << "Process row where first row is: " << desc[6] << std::endl;
                std::cout << "Process Col where first col is: " << desc[7] << std::endl;
                std::cout << "Leading Dimension: " << desc[8] << std::endl;
                std::cout << "Memory usage(data only) GB = "<<p.GBs<<std::endl;
        }
        return os;
}