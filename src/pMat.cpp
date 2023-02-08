#include "pMat.hpp"

using namespace ::std;

pMat::pMat()
{
    cout << "empty pMat created" << endl;
}
pMat::pMat(pMat *point)
{
    pG = point->pG;
    setupMat(point->M, point->N, point->type, point->block, point->cycles, 0.0, true);
}
pMat::pMat(int m, int n, PGrid *pGp)
{
    pG = pGp;
    setupMat(m, n, 0, 0, 1, 0.0, true);
}
pMat::pMat(int m, int n, PGrid *pGp, bool stdout)
{
    pG = pGp;
    setupMat(m, n, 0, 0, 1, 0.0, stdout);
}
pMat::pMat(int m, int n, PGrid *pGp, int t, int b, double init)
{
    pG = pGp;
    setupMat(m, n, t, b, 1, init, true);
}
pMat::pMat(int m, int n, PGrid *pGp, int t, int b, double init, bool stdout)
{
    pG = pGp;
    setupMat(m, n, t, b, 1, init, stdout);
}
pMat::pMat(int m, int n, PGrid *pGp, int t, int b, int c, double init)
{
    pG = pGp;
    setupMat(m, n, t, b, c, init, true);
}
pMat::pMat(int m, int n, PGrid *pGp, int t, int b, int c, double init, bool stdout)
{
    pG = pGp;
    setupMat(m, n, t, b, c, init, stdout);
}
pMat::~pMat()
{

}
void destroyPMat(pMat *A, bool stdout)
{
    if ((A->pG->rank == 0) && stdout)
        cout << "Deallocating Distributed Matrix" << endl;
    delete A;
}
void destroyPMat(pMat *A)
{
    destroyPMat(A, true);
}
void pMat::setupMat(int m, int n, int t, int b, int c, double init, bool stdout)
{
    M = m;
    N = n;
    type = t;
    block = b;
    cycles = c;
    printRank = pG->printRank;

    if (stdout)
    {
        cout << "Creating Matrix" << endl
            << "M = " << M << " N = " << N << endl;
    }
    if (block == 0) // square blocks
    {
        nb = 128;
        mb = nb;
        if (mb > (M / pG->getDim(1)))
            mb = max(1, M / (pG->getDim(1) * cycles));
        if (nb > (N / pG->getDim(0)))
            nb = max(1, N / (pG->getDim(0) * cycles));
        mb = min(mb, nb);
        nb = mb;
        if (stdout)
            cout << "mb/nb = " << nb << endl;
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else if (block == 1) // load blocks
    {
        nb = 128;
        mb = M;
        if (nb > (N / pG->getDim(0)))
            nb = max(1, N / pG->getDim(0));
        if (stdout)
        {
            cout << "mb is " << mb << endl;
            cout << "nb is " << nb << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else if (block == 2) // p0 block
    {
        nb = N;
        mb = M;
        if (stdout)
        {
            cout << "nb is " << nb << endl;
            cout << "mb is " << mb << endl;
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
    else if (block == 3) // synched block
    {
        nb = N / (pG->prow);
        mb = M / (pG->pcol);
    }

    myRC[0] = numroc_(&M, &mb, &(pG->myrow), &i_zero, &(pG->prow));
    myRC[1] = numroc_(&N, &nb, &(pG->mycol), &i_zero, &(pG->pcol));

    nelements = (long long)myRC[0] * (long long)myRC[1];
    for (int i = 0; i < 2; i++)
        myRC[i] = max(1, myRC[i]);

    if (type == 0)
    {
        if (stdout)
            cout << "Mat is Double" << endl;
        dataD.resize(nelements, init);
        MBs = nelements * 8.0 / (1.0e6);
    }
    else if (type == 1)
    {
        if (stdout)
            cout << "Mat is Complex" << endl;
        dataC.resize(nelements);
        for (int i = 0; i < nelements; i++)
        {
            dataC[i].real = init;
            dataC[i].imag = 0.0;
            MBs = nelements * 16.0 / (1.0e6);
        }
    }
    else
    {
        cout << "invalid type" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if (stdout)
    {
        cout << "nelements " << nelements << endl;
        cout << myRC[0] << " " << myRC[1] << endl;
    }

    desc = new int[9];
    descinit_(desc, &M, &N, &mb, &nb, &i_zero, &i_zero, &(pG->icntxt), &myRC[0], &info);
    if (info != 0)
    {
        cout << "Error in descriptor setup in argument, info=" << -info << endl;
    }

    if (stdout)
        cout << "Matrix Constructed" << endl;
}

void pMat::switchType(int t)
{
    if (type == t)
    {
        if (pG->printRank)
            cout << "matrix is already of data type " << type << endl;
    }
    else if (t == 0)
    {
        if (pG->printRank)
            cout << "Switching to Double " << endl;
        dataD.resize(nelements);
        for (int i = 0; i < nelements; i++)
            dataD[i] = dataC[i].real;
        dataC.clear();
        type = 0;
        if (pG->printRank)
            cout << "Matrix is now type" << type << endl;
    }
    else if (t == 1)
    {
        if (pG->printRank)
            cout << "Switching to Complex" << endl;

        dataC.resize(nelements);
        for (int i = 0; i < nelements; i++)
        {
            dataC[i].real = dataD[i];
            dataC[i].imag = 0;
        }
        dataD.clear();
        type = 1;
        if (pG->printRank)
            cout << "Matrix is now type " << type << endl;
    }
}

void pMat::printMat()
{
    for (int p = 0; p < pG->size; p++)
    {
        if (p == pG->rank)
        {
            cout << "processor " << p << " has :" << endl;
            for (int i = 0; i < nelements; i++)
            {
                if (i % myRC[0] == 0)
                    cout << endl;

                if (type == 0)
                    cout << i << " " << dataD[i] << endl;
                else if (type == 1)
                    cout << dataC[i].real << " + " << dataC[i].imag << " i " << endl;
            }
            fflush(stdout);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void pMat::write_bin(string filename)
{
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File fH = NULL;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fH);
    if (printRank)
    {
        MPI_File_write(fH, &M, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_write(fH, &N, 1, MPI_INT, MPI_STATUS_IGNORE);
        cout << "Write Start: " << filename << endl;
        cout << "M=" << M << ", mb=" << mb << ", N=" << N << ", nb=" << nb << endl;
    }
    MPI_File_close(&fH);
    MPI_Barrier(MPI_COMM_WORLD);
    int disp = 2 * sizeof(int);
    int dims[2] = {M, N};
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {mb, nb};
    MPI_Datatype darray;
    MPI_Type_create_darray(pG->size, pG->rank, 2, dims, distribs, dargs, pG->pdims, MPI_ORDER_FORTRAN, MPI_DOUBLE, &darray);
    MPI_Type_commit(&darray);
    int tsize, mpiEls;
    MPI_Type_size(darray, &tsize);
    mpiEls = tsize / (sizeof(double));
    if (nelements != mpiEls)
    {
        cout << "Allocation via MPI " << mpiEls << " and pblacs " << nelements << " is different" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    double t2, t1;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fH);
    if (printRank)
        cout << "MPI Allocation " << mpiEls << " , pblacs Allocation " << nelements << endl;
    t1 = MPI_Wtime();
    MPI_File_set_view(fH, disp, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_write_all(fH, dataD.data(), mpiEls, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fH);
    t2 = MPI_Wtime();
    if (printRank)
        cout << "Write time is " << t2 - t1 << endl;
}

void pMat::read_bin(string filename)
{
    int rM, rN;
    double t2, t1;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File fH;
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fH);
    MPI_File_read_all(fH, &rM, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_read_all(fH, &rN, 1, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fH);

    if ((rN != N) || (rM != M))
    {
        cout << "bin file and matrix do not have same dimension" << endl
            << "File M=" << rM << ", N=" << rN << endl
            << "Matrix M=" << M << ", N=" << N << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    int dims[2] = {M, N};
    cout << "M = " << M << " N = " << N << endl;
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {mb, nb};
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
        cout << "Allocation via MPI " << mpiEls << " and pblacs " << nelements << " is different" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    cout << "MPI Allocation " << mpiEls << " , pblacs Allocation " << nelements << endl;
    cout << "Read Starting" << endl;
    t1 = MPI_Wtime();
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fH);
    MPI_File_set_view(fH, disp, MPI_DOUBLE, darray, "native", MPI_INFO_NULL);
    MPI_File_read_all(fH, dataD.data(), mpiEls, MPI_DOUBLE, MPI_STATUS_IGNORE);
    MPI_File_close(&fH);
    t2 = MPI_Wtime();
    if (printRank)
        cout << "read of " << filename << " took " << t2 - t1 << " seconds" << endl;

}

bool pMat::check_bin_size(string filename, int &mM, int &mN)
{
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_File fH;
    if (!MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fH))
    {
        MPI_File_read_all(fH, &mM, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_read_all(fH, &mN, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fH);
        cout << filename << " has data with M: " << mM << " N: " << mN << endl;
        return true;
    }
    else
    {
        cout << filename << " not found" << endl;
        return false;
    }
}

void pMat::matrix_Product(char tA, char tB, int m, int n, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc)
{

    if ((A->type == 0) && (B->type == 0) && (type == 0))
    {
        int IA = ia + 1;
        int JA = ja + 1;
        int IB = ib + 1;
        int JB = jb + 1;
        int IC = ic + 1;
        int JC = jc + 1;
        pdgemm_(&tA, &tB, &m, &n, &k, &alpha, A->dataD.data(), &IA, &JA, A->desc, B->dataD.data(), &IB, &JB, B->desc, &beta, dataD.data(), &IC, &JC, desc);
    }
    else if ((A->type == 1) && (B->type == 1) && (type == 1))
    {
        int IA = ia + 1;
        int JA = ja + 1;
        int IB = ib + 1;
        int JB = jb + 1;
        int IC = ic + 1;
        int JC = jc + 1;
#ifdef USE_MKL
        MKL_Complex16 Ac, Bc;
#else
        complex16 Ac, Bc;
#endif
        Ac.real = alpha;
        Ac.imag = 0;
        Bc.real = beta;
        Bc.imag = 0;
        pzgemm_(&tA, &tB, &m, &n, &k, &Ac, A->dataC.data(), &IA, &JA, A->desc, B->dataC.data(), &IB, &JB, B->desc, &Bc, dataC.data(), &IC, &JC, desc);
    }
    else
    {
        if (printRank)
            cout << "Other matrix datatypes not supported yet" << endl;
    }
}

void pMat::ATA_lowMem(int m, int n, double alpha, pMat* A, int ia, int ja, int ic, int jc)
{
    int minMN = min(A->M - ia, A->N - ja);

    // if no memory limitation, just compute normal inner product
    if (minMN < MAX_OPERATION_SIZE)
    {
        this->matrix_Product('T', 'N', n - ja, m - ia, A->M, A, ia, ja, A, ia, ja, alpha, 0.0, ic, jc);
    }
    else
    {
        int numBlocks = ceil((float)minMN / (float)MAX_OPERATION_SIZE);
        int numRowsBlock = (A->M - ia) / numBlocks;
        int numColsBlock = (A->N - ja) / numBlocks;

        int numRows_AT, numCols_AT, numCols_A;
        int IA_AT, JA_AT, IA_A, JA_A;

        // row block index of A^T
        for (int i = 0; i < numBlocks; ++i)
        {
            IA_AT = ia + i * numRowsBlock;
            if (i == (numBlocks - 1))
                numRows_AT = A->M - i * numRowsBlock;
            else
                numRows_AT = numRowsBlock;

            // column block index of A^T
            for (int j = 0; j < numBlocks; ++j)
            {
                JA_AT = ja + j * numColsBlock;
                if (j == (numBlocks - 1))
                    numCols_AT = A->N - j * numColsBlock;
                else
                    numCols_AT = numColsBlock;

                // row block index of A
                for (int x = 0; x < numBlocks; ++x)
                {
                    IA_A = ia + x * numRowsBlock;

                    // column block index of A
                    for (int y = 0; y < numBlocks; ++y)
                    {
                        JA_A = ja + y * numColsBlock;
                        if (y == (numBlocks - 1))
                            numCols_A = A->N - y * numColsBlock;
                        else
                            numCols_A = numColsBlock;

                        this->matrix_Product('T', 'N', numCols_AT, numCols_A, numRows_AT, A, IA_AT, JA_AT, A, IA_A, JA_A, alpha, 1.0, JA_AT, JA_A);
                    }
                }
            }

        }

    }

}

// symmetric matrix product A^T*A or A*A^T
void pMat::matrix_Product_sym(char uplo, char trans, int n, int k, double alpha, pMat *A, int ia, int ja, double beta, int ic, int jc)
{

    int IA = ia + 1;
    int JA = ja + 1;
    int IC = ic + 1;
    int JC = jc + 1;
    pdsyrk_(&uplo, &trans, &n, &k, &alpha, A->dataD.data(), &IA, &JA, A->desc, &beta, dataD.data(), &IC, &JC, desc);

}

// matrix vector product A*X, where X is a subvector of pMat* B
void pMat::matrix_vec_product(char trans, int m, int n, double alpha, pMat *A, int ia, int ja, pMat *B, int ib, int jb,
                             double beta, int ic, int jc)
{

    int IA = ia + 1;
    int JA = ja + 1;
    int IB = ib + 1;
    int JB = jb + 1;
    int IC = ic + 1;
    int JC = jc + 1;
    int INCB = 1;
    int INCC = 1;
    pdgemv_(&trans, &m, &n, &alpha, A->dataD.data(), &IA, &JA, A->desc, B->dataD.data(), &IB, &JB, B->desc, &INCB, &beta, dataD.data(), &IC, &JC, desc, &INCC);

}

void pMat::matrix_Sum(char tA, int m, int n, pMat *A, int ia, int ja, int ib, int jb, double alpha, double beta)
{

    if ((A->type == 0) && (type == 0))
    {
        int IA = ia + 1;
        int JA = ja + 1;
        int IB = ib + 1;
        int JB = jb + 1;
        pdgeadd_(&tA, &m, &n, &alpha, A->dataD.data(), &IA, &JA, A->desc, &beta, dataD.data(), &IB, &JB, desc);
    }
    else if ((A->type == 1) && (type == 1))
    {
        int IA = ia + 1;
        int JA = ja + 1;
        int IB = ib + 1;
        int JB = jb + 1;
#ifdef USE_MKL
        MKL_Complex16 Ac, Bc;
#else
        complex16 Ac, Bc;
#endif
        Ac.real = alpha;
        Ac.imag = 0;
        Bc.real = beta;
        Bc.imag = 0;
        pzgeadd_(&tA, &m, &n, &Ac, A->dataC.data(), &IA, &JA, A->desc, &Bc, dataC.data(), &IB, &JB, desc);
    }
    else
    {
        cout << "Other Formats not supported yet" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

// Element-wise sub-matrix multiplication, scaled by alpha
void pMat::matrix_elem_mult(char tA, int m, int n, double alpha, pMat *A, int ia, int ja, int ic, int jc)
{

    if ((A->type == 0) && (type == 0))
    {
        int IA = ia + 1;
        int JA = ja + 1;
        int IC = ic + 1;
        int JC = jc + 1;
        pdgemul(tA, m, n, alpha, A->dataD.data(), IA, JA, A->desc, dataD.data(), IC, JC, desc);
    }
    else
    {
        cout << "Other formats not supported yet" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

// scales row or column of matrix by factor alpha
// idx indicates the row or column index (zero-indexed)
// if scaleRow == true, scales row, otherwise scales column
void pMat::scale_col_row(double alpha, int idx, bool scaleRow)
{
    int inc, len;
    int ix, jx;
    if (scaleRow)
    {
        len = N;
        inc = M;
        ix = idx + 1;
        jx = 1;
    }
    else
    {
        len = M;
        inc = 1;
        ix = 1;
        jx = idx + 1;
    }

    pdscal_(&len, &alpha, dataD.data(), &ix, &jx, desc, &inc);

}

void pMat::svd_run(pMat *U, pMat *VT, vector<double> &S)
{
    svd_run(M, N, 0, 0, U, VT, S, true);
}

void pMat::svd_run(int M, int N, int ia, int ja, pMat *U, pMat *VT, vector<double> &S)
{
    svd_run(M, N, ia, ja, U, VT, S, true);
}

void pMat::svd_run(int M, int N, int ia, int ja, pMat *U, pMat *VT, vector<double> &S, bool stdout)
{
    int info = 0;
    string computeFlag = "V";
    const char *JOBU = computeFlag.c_str(), *JOBVT = computeFlag.c_str();
    vector<double> WORK(1);
    int LWORK = -1;
    int IA = ia + 1;
    int JA = ja + 1;
    int i_one = 1;
    double t2, t1;
    pdgesvd_(JOBU, JOBVT, &M, &N, dataD.data(), &IA, &JA, desc, S.data(), U->dataD.data(), &IA, &JA, U->desc, VT->dataD.data(), &i_one, &i_one, VT->desc, WORK.data(), &LWORK, &info);
    if (stdout)
        cout << "WORK= " << WORK[0] << ", LWORK= " << LWORK << ", info= " << info << endl;

    if (info < 0)
    {
        cout << "Error in SVD setup in argument, info=" << -info << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
    LWORK = WORK[0];
    WORK.resize(LWORK);
    if (stdout)
        cout << "Work Allocated: " << LWORK / (1e6) * 8 << " MB per processor" << endl;

    //SVD run
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    pdgesvd_(JOBU, JOBVT, &M, &N, dataD.data(), &IA, &JA, desc, S.data(), U->dataD.data(), &IA, &JA, U->desc, VT->dataD.data(), &i_one, &i_one, VT->desc, WORK.data(), &LWORK, &info);
    t2 = MPI_Wtime();
    if (stdout)
        cout << "SVD complete in " << t2 - t1 << " seconds" << endl;
    WORK.resize(0);

}

void pMat::mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S, int modeStart, int modeEnd,
                  int mosStep, PGrid *procGrid)
{
    int IA = ia + 1;
    int JA = ja + 1;
    int i_one = 1;
    double t2, t1;
    cout << "starting MOS" << endl;
    int minMN = min(M, N);

    pMat *corMat;
    pMat *corMatp0;
    pMat *Vp0;
    pMat *V;

    if ((mosStep == 1) || (mosStep == 2))
    {
        corMat = new pMat(minMN, minMN, procGrid, 0, 0, 0.0);
        if (mosStep == 2)
        {
            corMatp0 = new pMat(minMN, minMN, procGrid, 0, 2, 0.0);
            Vp0 = new pMat(minMN, minMN, procGrid, 0, 2, 0.0);
        }
    }

    if (minMN == N)
    {
        // computing correlation matrix
        if (mosStep == 1)
        {
            cout << "Calculating correlation matrix" << endl;
            corMat->ATA_lowMem(minMN, minMN, 1.0, this, 0, 0, 0, 0);
            // corMat->matrix_Product('T', 'N', minMN, minMN, M, this, 0, 0, this, 0, 0, 1.0, 0.0, 0, 0);
            cout << "Correlation matrix calculated" << endl;

            corMat->write_bin("corMat.bin");
            delete corMat;
        }
        // computing singular values and right singular vectors
        else if (mosStep == 2)
        {
            // configure corMat on proc 0
            string corMatString = "corMat.bin";
            corMat->read_bin(corMatString);
            corMatp0->changeContext(corMat);
            delete corMat;

            cout << minMN << endl;
            if (procGrid->rank == 0)
            {
                cout << "Start eigensolve" << endl;
                t1 = MPI_Wtime();

                // map corMatp0 to Eigen data structure
                Eigen::MatrixXd corMatEig = Eigen::Map<Eigen::MatrixXd>(corMatp0->dataD.data(), minMN, minMN);

                // compute eigensolve
                Eigen::EigenSolver<Eigen::MatrixXd> es(corMatEig);

                Eigen::MatrixXcd VEig(minMN, minMN);
                Eigen::VectorXcd SEig(minMN);

                SEig = es.eigenvalues();
                VEig = es.eigenvectors();

                // map Eigen objects to double vectors
                Eigen::Map<Eigen::MatrixXd>(Vp0->dataD.data(), VEig.rows(), VEig.cols()) = VEig.real();
                Eigen::Map<Eigen::VectorXd>(S.data(), minMN) = SEig.real();

                cout << S.size() << endl;
                // singular values of A are square roots of eigenvalues of A^T*A
                for (int i = 0; i < S.size(); i++)
                {
                    if (S[i] < 0.0)
                        S[i] = 0.0;
                    S[i] = sqrt(S[i]);
                }
                t2 = MPI_Wtime();
                cout << "finish eigensolve in " << t2 - t1 << " seconds" << endl;

                // argsort singular values in descending order
                vector<size_t> idx(S.size());
                iota(idx.begin(), idx.end(), 0);
                stable_sort(idx.begin(), idx.end(), [&S](size_t i1, size_t i2) { return S[i1] < S[i2]; });
                reverse(idx.begin(), idx.end());

                // sort singular values and right singular vectors
                // TODO: this is not memory efficient, defnitely need to fix this
                vector<double> S_sort(minMN);
                vector<double> V_sort(minMN * minMN);
                for (int i = 0; i < minMN; ++i)
                {
                    S_sort[i] = S[idx[i]];
                    for (int j = 0; j < minMN; ++j)
                    {
                        V_sort[i * minMN + j] = Vp0->dataD[idx[i] * minMN + j];
                    }
                }

                for (int i = 0; i < minMN; ++i)
                {
                    S[i] = S_sort[i];
                    for (int j = 0; j < minMN; ++j)
                    {
                        Vp0->dataD[i * minMN + j] = V_sort[i * minMN + j];
                    }
                }

                // write singular values
                FILE *fid;
                fid = fopen("S.bin", "wb");
                fwrite(S.data(), sizeof(double), N, fid);
                fclose(fid);
            }
            MPI_Barrier(MPI_COMM_WORLD);

            // write right singular vectors
            V = new pMat(minMN, minMN, procGrid, 0, 0, 0.0);
            V->changeContext(Vp0);
            cout << "write V" << endl;
            V->write_bin("V.bin");

            delete corMatp0;
            delete Vp0;

        }
        else if (mosStep == 3)
        {
            // compute left singular vectors

            V = new pMat(VT->M, VT->N, VT->pG, 0, VT->block, 0.0);
            string Vstring = "V.bin";
            V->read_bin(Vstring);

            FILE *fid;
            fid = fopen("S.bin", "rb");
            size_t warn = fread(S.data(), sizeof(double), N, fid);

            t1 = MPI_Wtime();
            int modeCount = 0;
            for (int i = modeStart - 1; i < modeEnd; i++)
            {
                cout << "Processing left singular vector " << (i + 1) << endl;
                // U->matrix_Product('N','N',M,1,minMN,this,0,0,V,0,i,(1.0/S[i]),0.0,0,modeCount);
                U->matrix_vec_product('N', M, N, (1.0 / S[i]), this, 0, 0, V, 0, i, 0.0, 0, modeCount);
                modeCount++;
            }
            t2 = MPI_Wtime();
            cout << "MOS SVD complete in " << t2 - t1 << " seconds" << endl;

            VT->transpose(V);
            delete V;

        }
        else
        {
            cout << "Invalid value of mosStep: " << mosStep << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        cout << "min M not supported yet" << endl;
        MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

// mos_run for full run
void pMat::mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S)
{
    int modeStart = 1;
    int modeEnd = N;
    mos_run(M, N, ia, ja, U, VT, S, modeStart, modeEnd);
}

void pMat::mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S, int modeStart, int modeEnd)
{
    int IA = ia + 1;
    int JA = ja + 1;
    int i_one = 1;
    double t2, t1;
    cout << "starting MOS" << endl;
    int minMN = min(M, N);

    pMat *corMat = new pMat(minMN, minMN, pG, 0, 0, 0.0, false);
    pMat *corMatp0 = new pMat(minMN, minMN, pG, 0, 2, 0.0, false);
    pMat *Vp0 = new pMat(minMN, minMN, pG, 0, 2, 0.0, false);

    if (minMN == N)
    {
        corMat->ATA_lowMem(minMN, minMN, 1.0, this, 0, 0, 0, 0);
        // corMat->matrix_Product('T', 'N', minMN, minMN, M, this, 0, 0, this, 0, 0, 1.0, 0.0, 0, 0);

        corMatp0->changeContext(corMat, false);
        delete corMat;

        // cout << "start eigensolve" << endl;
        if (pG->rank == 0)
        {
            // t1 = MPI_Wtime();
            Eigen::MatrixXd corMatEig = Eigen::Map<Eigen::MatrixXd>(corMatp0->dataD.data(), minMN, minMN);
            Eigen::EigenSolver<Eigen::MatrixXd> es(corMatEig);

            Eigen::MatrixXcd VEig(minMN, minMN);
            Eigen::VectorXcd SEig(minMN);

            // eigen gives eigenvalues in ascending order, reverse for descending order
            SEig = es.eigenvalues();
            VEig = es.eigenvectors();

            // map Eigen objects to double vectors
            Eigen::Map<Eigen::MatrixXd>(Vp0->dataD.data(), VEig.rows(), VEig.cols()) = VEig.real();
            Eigen::Map<Eigen::VectorXd>(S.data(), minMN) = SEig.real();

            // singular values of A are square roots of eigenvalues of A^T*A
            for (int i = 0; i < S.size(); i++)
            {
                if (S[i] < 0.0)
                        S[i] = 0.0;
                S[i] = sqrt(S[i]);
            }

            // argsort singular values in ascending order
            vector<size_t> idx(S.size());
            iota(idx.begin(), idx.end(), 0);
            stable_sort(idx.begin(), idx.end(), [&S](size_t i1, size_t i2) { return S[i1] < S[i2]; });

            reverse(idx.begin(), idx.end()); // reverse order to get singular values in descending order

            // sort singular values and right singular vectors
            // TODO: this is not memory efficient, defnitely need to fix this
            vector<double> S_sort(minMN);
            vector<double> V_sort(minMN * minMN);
            for (int i = 0; i < minMN; ++i)
            {
                S_sort[i] = S[idx[i]];
                for (int j = 0; j < minMN; ++j)
                {
                    V_sort[i * minMN + j] = Vp0->dataD[idx[i] * minMN + j];
                }
            }

            for (int i = 0; i < minMN; ++i)
            {
                S[i] = S_sort[i];
                for (int j = 0; j < minMN; ++j)
                {
                    Vp0->dataD[i * minMN + j] = V_sort[i * minMN + j];
                }
            }
        }

        delete corMatp0;

        // map V from rank 0 to distributed matrix
        MPI_Bcast(S.data(), S.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
        pMat *V = new pMat(VT->N, VT->M, VT->pG, 0, VT->block, 0.0,false);
        V->changeContext(Vp0, false);
        delete Vp0;

        VT->transpose(V);
        // delete V;

        int modeCount = 0;
        for (int i = modeStart - 1; i < modeEnd; i++)
        {
            // cout << modeCount
            // U->matrix_Product('N', 'T', M, 1, minMN, this, 0, 0, VT, i, 0, (1.0 / S[i]), 0.0, 0, modeCount);
            U->matrix_vec_product('N', M, N, (1.0 / S[i]), this, 0, 0, V, 0, i, 0.0, 0, modeCount);
            modeCount++;
        }

    }
    else
    {
            cout << "min M not supported yet" << endl;
            MPI_Abort(MPI_COMM_WORLD, -1);
    }
}

void pMat::qr_run(int m, int n, int ia, int ja, vector<int> &ipiv)
{
    qr_run(m, n, ia, ja, ipiv, "./", "P.bin", true);
}

void pMat::qr_run(int m, int n, int ia, int ja, vector<int> &ipiv, string outdir, bool stdout) {
	qr_run(m, n, ia, ja, ipiv, outdir, "P.bin", stdout);
}

void pMat::qr_run(int m, int n, int ia, int ja, vector<int> &ipiv, string outdir, string outfile, bool stdout)
{
    if (stdout)
        cout << "QR initializing" << endl;
    int IA = ia + 1;
    int JA = ja + 1;

    int JAnpm1 = JA + n - 1;
    int JAminMNm1 = JA + min(m, n) - 1;

    int ipiv_LOCc = max(1, numroc_(&JAnpm1, &nb, &(pG->mycol), &(desc[7]), &(pG->pcol)));
    int tau_LOCc = max(1, numroc_(&JAminMNm1, &nb, &(pG->mycol), &(desc[7]), &(pG->pcol)));
    if (stdout)
    {
        cout << "ipiv_LOCc= " << ipiv_LOCc << endl;
        cout << "tau_LOCc= " << tau_LOCc << endl;
    }
    vector<double> tau(1);
    vector<double> WORK(1);

    ipiv.resize(ipiv_LOCc);
    tau.resize(tau_LOCc);

    double t1, t2;
    int info = 0;
    int LWORK = -1;

    pdgeqpf_(&m, &n, dataD.data(), &IA, &JA, desc, ipiv.data(), tau.data(), WORK.data(), &LWORK, &info);
    if (stdout)
        cout << "WORK=" << WORK[0] << " ,LWORK= " << LWORK << ", info= " << info << endl;
    LWORK = WORK[0];
    WORK.resize(LWORK);
    if (stdout)
        cout << "WORK allocated starting QR" << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    t1 = MPI_Wtime();
    pdgeqpf_(&m, &n, dataD.data(), &IA, &JA, desc, ipiv.data(), tau.data(), WORK.data(), &LWORK, &info);
    MPI_Barrier(MPI_COMM_WORLD);
    t2 = MPI_Wtime();
    if (stdout)
        cout << "QR complete in " << t2 - t1 << " seconds info =" << info << endl;
    MPI_Barrier(MPI_COMM_WORLD);
    tau.resize(0);
    WORK.resize(0);

    MPI_Barrier(MPI_COMM_WORLD);
    int ONE = 1;
    MPI_File fH;
    string pivot_name = outdir + "/" + outfile;
    MPI_File_open(MPI_COMM_WORLD, pivot_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fH);
    if (printRank)
    {
        MPI_File_write(fH, &ONE, 1, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_write(fH, &N, 1, MPI_INT, MPI_STATUS_IGNORE);
        if (stdout)
        {
            cout << "Write Start" << endl;
            cout << "M=" << ONE << ", mb=" << mb << ", N=" << N << ", nb=" << nb << endl;
        }
    }
    MPI_File_close(&fH);
    MPI_Barrier(MPI_COMM_WORLD);
    int disp = 2 * sizeof(int);
    int dims[2] = {ONE, N};
    int distribs[2] = {MPI_DISTRIBUTE_CYCLIC, MPI_DISTRIBUTE_CYCLIC};
    int dargs[2] = {mb, nb};
    MPI_Datatype darray;
    MPI_Type_create_darray(pG->size, pG->rank, 2, dims, distribs, dargs, pG->pdims, MPI_ORDER_FORTRAN, MPI_INT, &darray);
    MPI_Type_commit(&darray);
    int tsize, mpiEls;
    MPI_Type_size(darray, &tsize);
    mpiEls = tsize / (sizeof(int));
    if (myRC[1] != mpiEls)
    {
        cout << "Allocation via MPI " << mpiEls << " and pblacs " << myRC[1] << " is different" << endl;
    }
    MPI_File_open(MPI_COMM_WORLD, pivot_name.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fH);
    if (printRank && stdout)
            cout << "MPI Allocation " << mpiEls << " , pblacs Allocation " << myRC[1] << endl;
    t1 = MPI_Wtime();
    MPI_File_set_view(fH, disp, MPI_INT, darray, "native", MPI_INFO_NULL);
    MPI_File_write_all(fH, ipiv.data(), mpiEls, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fH);
    t2 = MPI_Wtime();
    if (printRank && stdout)
        cout << "Write time is " << t2 - t1 << endl;

}

void pMat::transpose(pMat *A)
{
    transpose(A, this->M, this->N, 0, 0);
}

void pMat::transpose(pMat *A, int m, int n, int ia, int ja)
{
    int IA = ia + 1;
    int JA = ja + 1;
    //cout << "Copying transpose" << endl;

    if ((m != M) && (n != A->M))
    {
        cout << "transpose dimension mismatch" << endl;
        cout << m << " " << M << " " << n << " " << A->M << endl;
        throw(-1);
    }
    int i_one = 1;
    double ONE = 1.0;
    double ZERO = 0.0;
    pdtran_(&m, &n, &ONE, A->dataD.data(), &IA, &JA, A->desc, &ZERO, dataD.data(), &i_one, &i_one, desc);
}

void pMat::changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb, bool stdout)
{
    int IA = ia + 1;
    int JA = ja + 1;
    int IB = ib + 1;
    int JB = jb + 1;
    if (stdout)
        cout << "Copying Matrix" << endl
            << "m = " << m << " , "
            << "n = " << n << endl;
    int i_one = 1;
    if (type == 0)
    {
        if (stdout)
            cout << "Double changed pGrid" << endl;
        pdgemr2d_(&m, &n, A->dataD.data(), &IA, &JA, A->desc, dataD.data(), &IB, &JB, desc, &(pG->icntxt));
    }
    if (type == 1)
    {
        if (stdout)
            cout << "Complex changed pGrid" << endl;
        pzgemr2d_(&m, &n, A->dataC.data(), &IA, &JA, A->desc, dataC.data(), &IB, &JB, desc, &(pG->icntxt));
    }
}

void pMat::changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb)
{
    changeContext(A, m, n, ia, ja, ib, jb, true);
}

void pMat::changeContext(pMat *A)
{
    changeContext(A, M, N, 0, 0, 0, 0, true);
}

void pMat::changeContext(pMat *A, bool stdout)
{
    changeContext(A, M, N, 0, 0, 0, 0, stdout);
}

void pMat::dMax(int dim, int rc, double &val, int &index)
{

    if (dim == 0)
    {
        int IA = 1, JA = rc + 1, i_one = 1;
        pdamax_(&M, &val, &index, dataD.data(), &IA, &JA, desc, &i_one);
    }
    if (dim == 1)
    {
        int IA = rc + 1, JA = 1, i_one = 1;
        pdamax_(&N, &val, &index, dataD.data(), &IA, &JA, desc, &M);
    }

    index--; // return to C indexing

}

// Computes global argmax of one-dimensional pMat
int pMat::argmax_vec()
{
    // only valid for "vectors" (one-dimensiona pMat)
    int vecType;
    if (N == 1) {
        vecType = 0; // column vector
    } else if (M == 1) {
        vecType = 1; // row vector
    } else {
        cout << "pMat::argmax_vec only accepts one-dimensiona pMats" << endl;
        cout << "This pMat has dimensions (" << M << ", " << N << ")" << endl;
        throw(-1);
    }

    double maxLocal = 0.0;
    double maxGlobal = 0.0;
    int argMaxLocal = -1;
    int argMaxGlobal = -1;

    this->dMax(vecType, 0, maxLocal, argMaxLocal);

    // Allreduce(MPI_MAX) the maximum value
    MPI_Allreduce(&maxLocal, &maxGlobal, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // if rank doesn't own the max, set argMaxLocal = -1 so that MPI_MAX can determine the correct argMaxGlobal
    if (maxLocal != maxGlobal) {
        argMaxLocal = -1;
    }
    MPI_Allreduce(&argMaxLocal, &argMaxGlobal, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    return argMaxGlobal;
}

void pMat::dSum(int dim, int rc, double &val)
{

    if (dim == 0)
    {
        int IA = 1, JA = rc + 1, i_one = 1;
        pdasum_(&M, &val, dataD.data(), &IA, &JA, desc, &i_one);
    }
    if (dim == 1)
    {
        int IA = rc + 1, JA = 1, i_one = 1;
        pdasum_(&N, &val, dataD.data(), &IA, &JA, desc, &M);
    }
}

void pMat::pinv(pMat *A)
{
    pMat *UU, *VV;
    UU = new pMat(A->M, min(A->M, A->N), A->pG, false);
    VV = new pMat(min(A->M, A->N), A->N, A->pG, false);
    vector<double> SS(min(A->M, A->N), 0.0);

    A->svd_run(A->M, A->N, 0, 0, UU, VV, SS, false);
    this->matrix_Product('T', 'T', VV->N, UU->M, 1, VV, 0, 0, UU, 0, 0, 1.0 / SS[0], 0.0, 0, 0);
    for (int i = 1; i < SS.size(); i++)
    {
        if (SS[i] > numeric_limits<double>::epsilon() * max(A->M, A->N) * SS[0])
            this->matrix_Product('T', 'T', VV->N, UU->M, 1, VV, i, 0, UU, 0, i, 1.0 / SS[i], 1.0, 0, 0);
    }
    delete UU;
    delete VV;
}

// solve over-/under-determined linear system AX = B
// on exit, solutions are written to columns of B
// on exit, A is overwritten with QR decomposition info (pretty much destroyed)

void pMat::leastSquares(pMat *A)
{
    leastSquares('N', M, N, A->N, A, 0, 0, 0, 0);
}

void pMat::leastSquares(char trans, int m, int n, int nrhs, pMat *A, int ia, int ja, int ib, int jb)
{

	int info = 0;
	vector<double> WORK(1);
	int LWORK = -1;
	int IA = ia + 1;
	int JA = ja + 1;
	int IB = ib + 1;
	int JB = jb + 1;

	// get LWORK and WORK
	pdgels_(&trans, &m, &n, &nrhs, A->dataD.data(), &IA, &JA, A->desc, dataD.data(), &IB, &JB, desc, WORK.data(), &LWORK, &info);

	if (info < 0) {
		cout << "Error in least-squares solve setup in argument: " << -info << endl;
		throw(-1);
	}

	// set up real run
	LWORK = WORK[0];
	WORK.resize(LWORK);

	// least squares solve
	double t1, t2;
	t1 = MPI_Wtime();
	pdgels_(&trans, &m, &n, &nrhs, A->dataD.data(), &IA, &JA, A->desc, dataD.data(), &IB, &JB, desc, WORK.data(), &LWORK, &info);
	t2 = MPI_Wtime();
	WORK.resize(0);

}

void pMat::outerProductSum(pMat *U, char UT, pMat *VT, char VTT, vector<double> &S, int inv)
{
    if (inv == 1)
    {
        this->matrix_Product(UT, VTT, U->M, VT->N, 1, U, 0, 0, VT, 0, 0, 1 / S[0], 0.0, 0, 0);
        for (int i = 1; i < S.size(); i++)
        {
            this->matrix_Product(UT, VTT, U->M, VT->N, 1, U, 0, i, VT, i, 0, 1 / S[i], 1.0, 0, 0);
        }
    }
    else
    {
        this->matrix_Product(UT, VTT, U->M, VT->N, 1, U, 0, 0, VT, 0, 0, S[0], 0.0, 0, 0);
        for (int i = 1; i < S.size(); i++)
        {
            this->matrix_Product(UT, VTT, U->M, VT->N, 1, U, 0, i, VT, i, 0, S[i], 1.0, 0, 0);
        }
    }
}

ostream &operator<<(ostream &os, const pMat &p)
{
    cout << "Descriptor type: " << p.desc[0] << endl;
    cout << "BLACS context: " << p.desc[1] << endl;
    cout << "Global Rows: " << p.desc[2] << endl;
    cout << "Global Cols: " << p.desc[3] << endl;
    cout << "Row Blocking factor: " << p.desc[4] << endl;
    cout << "Column Blocking factor: " << p.desc[5] << endl;
    cout << "Process row where first row is: " << p.desc[6] << endl;
    cout << "Process Col where first col is: " << p.desc[7] << endl;
    cout << "Leading Dimension: " << p.desc[8] << endl;
    cout << "Memory usage(data only) MB = " << p.MBs << endl;
    return os;
}

// determine index of element in dataD corresponding to global indices (I, J)
// returns -1 if process does not own global index
// zero-indexed
int pMat::getDataIndex(int I, int J)
{
    int l = I / (pG->prow * mb);
    int m = J / (pG->pcol * nb);

    int x = I % mb;
    int y = J % nb;

    // TODO: should this be int64?
    int localIndex = -1;
    if ((pG->myrow == (I / mb) % pG->prow) && (pG->mycol == (J / nb) % pG->pcol))
    {
        localIndex = (m * nb + y) * myRC[0] + l * mb + x;
    }

    return localIndex;
}

// retrieve element from pMat, no matter what process owns the element
double pMat::getElement(int I, int J)
{
    int index = getDataIndex(I, J);
    double temp = 0.0;
    if (index >= 0)
    {
        assert(index < nelements);
        temp = dataD[index];
    }
    MPI_Allreduce(MPI_IN_PLACE, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return temp;
}

// retrieve element from pMat, only if the process owns the element
// otherwise returns input value of temp (zero by default)
double pMat::getLocalElement(int I, int J)
{
    return getLocalElement(I, J, 0.0);
}

double pMat::getLocalElement(int I, int J, double temp)
{
    int index = getDataIndex(I, J);
    if (index >= 0)
    {
        assert(index < nelements);
        temp = dataD[index];
    }
    return temp;
}

void pMat::setElement(int I, int J, double val)
{
    setElement(I, J, val, true);
}

void pMat::setElement(int I, int J, double val, bool barrier)
{
        int index = getDataIndex(I, J);
        if (index >= 0)
        {
            assert(index < nelements);
            dataD[index] = val;
        }
        if (barrier)
        {
            MPI_Barrier(MPI_COMM_WORLD);
        }
}

bool operator==(pMat const &p1, pMat const &p2)
{
    if ((p1.M != p2.M) || (p1.N != p2.N))
    {
        cout << "Dim mismatch" << endl;
        return false;
    }
    else if ((p1.nelements != p2.nelements))
    {
        cout << "element mismatch" << endl;
        return false;
    }
    else
    {
        if (p1.type == 0)
        {
            for (int i = 0; i < p1.nelements; i++)
            {
                double epsilon = .01;
                if (p1.dataD[i] > 0)
                {
                    if ((p1.dataD[i] + p1.dataD[i] * epsilon) < p2.dataD[i] || (p1.dataD[i] - p1.dataD[i] * epsilon) > p2.dataD[i])
                    {
                        cout << "element " << i << " does not match" << endl;
                        cout << (p1.dataD[i] - p1.dataD[i] * epsilon) << "< " << p2.dataD[i] << " < " << (p1.dataD[i] + p1.dataD[i] * epsilon) << endl;
                        cout << p1.dataD[i] << "!=" << p2.dataD[i] << endl;
                        return false;
                    }
                }
                else
                {
                    if ((p1.dataD[i] + p1.dataD[i] * epsilon) > p2.dataD[i] || (p1.dataD[i] - p1.dataD[i] * epsilon) < p2.dataD[i])
                    {
                        cout << "element " << i << " does not match" << endl;
                        cout << (p1.dataD[i] - p1.dataD[i] * epsilon) << "< " << p2.dataD[i] << " < " << (p1.dataD[i] + p1.dataD[i] * epsilon) << endl;
                        cout << p1.dataD[i] << "!=" << p2.dataD[i] << endl;
                        return false;
                    }
                }
            }
        }
        else if (p1.type == 1)
        {
            for (int i = 0; i < p1.nelements; i++)
            {
                if ((p1.dataC[i].real != p2.dataC[i].real) || (p1.dataC[i].imag != p2.dataC[i].imag))
                    return false;
            }
        }
    }

    return true;
}

void pMat::commCreate(MPI_Comm &col_comm, int dim)
{
    if (dim == 0)
    {
        MPI_Comm_split(MPI_COMM_WORLD, pG->mycol, pG->rank, &col_comm);
    }
    else if (dim == 1)
    {
        MPI_Comm_split(MPI_COMM_WORLD, pG->myrow, pG->rank, &col_comm);
    }
    else
        MPI_Abort(MPI_COMM_WORLD, -1);
}
