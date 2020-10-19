#include "pMat.hpp"

using namespace ::std;

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
        if (type == 0)
                MPI_Dims_create(size, 2, pdims);
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
        cout << rank << " " << size;
        cout << "initializing Cblacs" << endl;
        Cblacs_pinfo(&rank, &size);
        Cblacs_get(-1, 0, &icntxt);
        Cblacs_gridinit(&icntxt, order, pdims[0], pdims[1]);
        Cblacs_gridinfo(icntxt, &pcol, &prow, &myrow, &mycol);
        cout << "Processor Grid M(cols)=" << pdims[0] << " N(rows)=" << pdims[1] << endl;
        cout << "local processor row" << myrow << " and column " << mycol << endl;
        cout << "PBLACS is row major by default and everything is else is Column major (because reasons)" << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        delete[] order;
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
pMat::pMat()
{
        cout << "empty pMat created" << endl;
}
pMat::pMat(pMat *point)
{
        pG = point->pG;
        setupMat(point->M, point->N, point->type, point->block, point->cycles, 0.0);
}
pMat::pMat(int m, int n, PGrid *pGp)
{
        pG = pGp;
        setupMat(m, n, 0, 0, 1, 0.0);
}
pMat::pMat(int m, int n, PGrid *pGp, int t, int b, double init)
{
        pG = pGp;
        setupMat(m, n, t, b, 1, init);
}
pMat::pMat(int m, int n, PGrid *pGp, int t, int b, int c, double init)
{
        pG = pGp;
        setupMat(m, n, t, b, c, init);
}
pMat::~pMat()
{
        if (pG->rank == 0)
                cout << "Deallocating Distributed Matrix" << endl;
}
void pMat::setupMat(int m, int n, int t, int b, int c, double init)
{
        M = m;
        N = n;
        type = t;
        block = b;
        cycles = c;
        printRank = pG->printRank;

        cout << "Creating Matrix" << endl
             << "M=" << M << " N=" << N << endl;
        if (block == 0) //square blocks
        {
                nb = 128;
                mb = nb;
                if (mb > (M / pG->getDim(1)))
                        mb = std::max(1, M / (pG->getDim(1) * cycles));
                if (nb > (N / pG->getDim(0)))
                        nb = std::max(1, N / (pG->getDim(0) * cycles));
                mb = std::min(mb, nb);
                nb = mb;
                cout << "mb/nb = " << nb << endl;
                MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (block == 1) //load blocks
        {
                nb = 128;
                mb = M;
                if (nb > (N / pG->getDim(0)))
                        nb = std::max(1, N / pG->getDim(0));
                cout << "mb is " << mb << endl;
                cout << "nb is " << nb << endl;
                MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (block == 2) //p0 block
        {
                nb = N;
                mb = M;
                cout << "nb is " << nb << endl;
                cout << "mb is " << mb << endl;
                MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (block == 3) //synched block
        {
                nb = N / (pG->prow);
                mb = M / (pG->pcol);
        }
        myRC[0] = numroc_(&M, &mb, &(pG->myrow), &i_zero, &(pG->prow));
        //cout<<myRC[0]<<" "<<M<<" "<<mb<<" "<<pG->myrow<<" "<<pG->prow<<endl;

        myRC[1] = numroc_(&N, &nb, &(pG->mycol), &i_zero, &(pG->pcol));
        //cout<<myRC[1]<<" "<<N<<" "<<nb<<" "<<pG->mycol<<" "<<pG->pcol<<endl;

        //if (myRC[0] < 1)
        //        myRC[0] = 1;
        //if (myRC[1] < 1)
        //        myRC[1] = 1;
        nelements = (long long)myRC[0] * (long long)myRC[1];
        for (int i = 0; i < 2; i++)
                myRC[i] = std::max(1, myRC[i]);

        if (type == 0)
        {
                cout << "Mat is Double" << endl;
                dataD.resize(nelements, init);
                MBs = nelements * 8.0 / (1.0e6);
        }
        else if (type == 1)
        {
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
                throw(-1);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        cout << "nelements " << nelements << endl;
        cout << myRC[0] << " " << myRC[1] << endl;

        descinit_(desc, &M, &N, &mb, &nb, &i_zero, &i_zero, &(pG->icntxt), &myRC[0], &info);
        //for(int i=0;i<9;i++)
        //        cout<<desc[i]<<endl;
        if (info != 0)
        {
                cout << "Error in descriptor setup in argument, info=" << -info << endl;
        }
        cout << "Matrix Constructed" << endl;
}

void pMat::switchType(int t)
{
        if (type == t)
        {
                if (pG->printRank)
                        cout << "matrix is already of data type " << type << endl;
                ;
                return;
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
                ;
                return;
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
                ;
                return;
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

int pMat::write_bin(std::string filename)
{
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_File fH;
        MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fH);
        if (printRank)
        {
                MPI_File_write(fH, &M, 1, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_write(fH, &N, 1, MPI_INT, MPI_STATUS_IGNORE);
                cout << "Write Start: " << filename << endl;
                cout << "M=" << M << "mb=" << mb << "N=" << N << "nb=" << nb << endl;
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
                throw(-1);
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

int pMat::read_bin(string &filename)
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
                throw(-1);
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
                throw(-1);
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

        return 1;
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

int pMat::matrix_Product(char tA, char tB, int m, int n, int k, pMat *A, int ia, int ja, pMat *B, int ib, int jb, double alpha, double beta, int ic, int jc)
{

        if ((A->type == 0) && (B->type == 0) && (type == 0))
        {
                int IA = ia + 1;
                int JA = ja + 1;
                int IB = ib + 1;
                int JB = jb + 1;
                int IC = ic + 1;
                int JC = jc + 1;
                pdgemm(&tA, &tB, &m, &n, &k, &alpha, A->dataD.data(), &IA, &JA, A->desc, B->dataD.data(), &IB, &JB, B->desc, &beta, dataD.data(), &IC, &JC, desc);
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
                pzgemm(&tA, &tB, &m, &n, &k, &Ac, A->dataC.data(), &IA, &JA, A->desc, B->dataC.data(), &IB, &JB, B->desc, &Bc, dataC.data(), &IC, &JC, desc);
        }
        else
        {
                if (printRank)
                        cout << "Other matrix datatypes not supported yet" << endl;
        }
        return 0;
}

// symmetric matrix product A^T*A or A*A^T
int pMat::matrix_Product_sym(char uplo, char trans, int n, int k, double alpha, pMat *A, int ia, int ja, double beta, int ic, int jc)
{

        int IA = ia + 1;
        int JA = ja + 1;
        int IC = ic + 1;
        int JC = jc + 1;
        pdsyrk(&uplo, &trans, &n, &k, &alpha, A->dataD.data(), &IA, &JA, A->desc, &beta, dataD.data(), &IC, &JC, desc);

        return 0;
}

// matrix vector product A*X, where X is a subvector of pMat* B
int pMat::matrix_vec_product(char trans, int m, int n, double alpha, pMat *A, int ia, int ja, pMat *B, int ib, int jb,
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
        pdgemv(&trans, &m, &n, &alpha, A->dataD.data(), &IA, &JA, A->desc, B->dataD.data(), &IB, &JB, B->desc, &INCB, &beta, dataD.data(), &IC, &JC, desc, &INCC);

        return 0;
}

int pMat::matrix_Sum(char tA, int m, int n, pMat *A, int ia, int ja, int ib, int jb, double alpha, double beta)
{

        if ((A->type == 0) && (type == 0))
        {
                if (printRank)
                        cout << "Double Sum" << endl;
                int IA = ia + 1;
                int JA = ja + 1;
                int IB = ib + 1;
                int JB = jb + 1;
                pdgeadd(&tA, &m, &n, &alpha, A->dataD.data(), &IA, &JA, A->desc, &beta, dataD.data(), &IB, &JB, desc);
        }
        else if ((A->type == 1) && (type == 1))
        {
                if (printRank)
                        cout << "Complex Sum" << endl;
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
                pzgeadd(&tA, &m, &n, &Ac, A->dataC.data(), &IA, &JA, A->desc, &Bc, dataC.data(), &IB, &JB, desc);
        }
        else
        {
                if (printRank)
                        cout << "Other Formats not supported yet" << endl;
        }
        return 0;
}

int pMat::svd_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S)
{
        int info = 0;
        string computeFlag = "V";
        const char *JOBU = computeFlag.c_str(), *JOBVT = computeFlag.c_str();
        std::vector<double> WORK(1);
        int LWORK = -1;
        int IA = ia + 1;
        int JA = ja + 1;
        int i_one = 1;
        double t2, t1;
        pdgesvd(JOBU, JOBVT, &M, &N, dataD.data(), &IA, &JA, desc, S.data(), U->dataD.data(), &IA, &JA, U->desc, VT->dataD.data(), &i_one, &i_one, VT->desc, WORK.data(), &LWORK, &info);

        cout << "WORK= " << WORK[0] << ", LWORK= " << LWORK << ", info= " << info << endl;

        if (info < 0)
        {
                cout << "Error in SVD setup in argument, info=" << -info << endl;
                throw(-1);
        }
        LWORK = WORK[0];
        WORK.resize(LWORK);
        cout << "Work Allocated: " << LWORK / (1e6) * 8 << " MB per processor" << endl;
        //SVD run
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        pdgesvd(JOBU, JOBVT, &M, &N, dataD.data(), &IA, &JA, desc, S.data(), U->dataD.data(), &IA, &JA, U->desc, VT->dataD.data(), &i_one, &i_one, VT->desc, WORK.data(), &LWORK, &info);
        t2 = MPI_Wtime();
        cout << "SVD complete in " << t2 - t1 << " seconds" << endl;
        WORK.resize(0);

        return 1;
}

int pMat::mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S, int modeStart, int modeEnd,
                  int mosStep, PGrid *procGrid)
{
        int IA = ia + 1;
        int JA = ja + 1;
        int i_one = 1;
        double t2, t1;
        cout << "starting MOS" << endl;
        int minMN = std::min(M, N);

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
                        // corMat->matrix_Product_sym('U', 'T', minMN, M, 1.0, this, 0, 0, 0.0, 0, 0);
                        corMat->matrix_Product('T', 'N', minMN, minMN, M, this, 0, 0, this, 0, 0, 1.0, 0.0, 0, 0);
                        cout << "Correlation matrix calculated" << endl;

                        corMat->write_bin("corMat.bin");
                        delete corMat;
                        return 0;
                }

                // computing singular values and right singular vectors
                if (mosStep == 2)
                {

                        // configure corMat on proc 0
                        std::string corMatString = "corMat.bin";
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
                                // corMatEig = corMatEig.selfadjointView<Eigen::Upper>();

                                // compute eigensolve
                                // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(corMatEig);
                                Eigen::EigenSolver<Eigen::MatrixXd> es(corMatEig);

                                Eigen::MatrixXcd VEig(minMN, minMN);
                                Eigen::VectorXcd SEig(minMN);
                                // Eigen::MatrixXd VEig(minMN,minMN);
                                // Eigen::VectorXd SEig(minMN);

                                // Eigen gives eigenvalues in ascending order, reverse for descending order
                                // SEig = es.eigenvalues().reverse();
                                // VEig = es.eigenvectors().rowwise().reverse();
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
                                        S[i] = std::sqrt(S[i]);
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

                        return 0;
                }
                else if (mosStep == 3)
                {
                        // compute left singular vectors

                        V = new pMat(VT->M, VT->N, VT->pG, 0, VT->block, 0.0);
                        std::string Vstring = "V.bin";
                        V->read_bin(Vstring);

                        FILE *fid;
                        fid = fopen("S.bin", "rb");
                        fread(S.data(), sizeof(double), N, fid);

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

                        return 0;
                }

                cout << "Invalid value of mosStep: " << mosStep << endl;
                throw(-1);
        }

        else
        {
                cout << "min M not supported yet" << endl;
                throw(-1);
        }
}

// mos_run for full run
int pMat::mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S)
{

        int modeStart = 1;
        int modeEnd = N;
        mos_run(M, N, ia, ja, U, VT, S, modeStart, modeEnd);
}

int pMat::mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S, int modeStart, int modeEnd)
{
        int IA = ia + 1;
        int JA = ja + 1;
        int i_one = 1;
        double t2, t1;
        cout << "starting MOS" << endl;
        int minMN = std::min(M, N);

        pMat *corMat = new pMat(minMN, minMN, pG, 0, 0, 0.0);
        pMat *corMatp0 = new pMat(minMN, minMN, pG, 0, 2, 0.0);
        pMat *Vp0 = new pMat(minMN, minMN, pG, 0, 2, 0.0);

        if (minMN == N)
        {
                corMat->matrix_Product('T', 'N', minMN, minMN, M, this, 0, 0, this, 0, 0, 1.0, 0.0, 0, 0);
                cout << "cor Mat created" << endl;

                corMatp0->changeContext(corMat);
                delete corMat;

                cout << "start eigensolve" << endl;
                if (pG->rank == 0)
                {
                        t1 = MPI_Wtime();
                        Eigen::MatrixXd corMatEig = Eigen::Map<Eigen::MatrixXd>(corMatp0->dataD.data(), minMN, minMN);
                        // corMatEig = corMatEig.selfadjointView<Eigen::Upper>();
                        // Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(corMatEig);
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
                                S[i] = std::sqrt(S[i]);
                        }
                        t2 = MPI_Wtime();
                        cout << "Finish eigensolve in " << t2 - t1 << " seconds" << endl;

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
                pMat *V = new pMat(VT->M, VT->N, VT->pG, 0, VT->block, 0.0);
                V->changeContext(Vp0);
                delete Vp0;

                VT->transpose(V);
                delete V;

                t1 = MPI_Wtime();
                int modeCount = 0;
                for (int i = modeStart - 1; i < modeEnd; i++)
                {
                        cout << "Processing left singular vector " << (i + 1) << endl;
                        U->matrix_Product('N', 'T', M, 1, minMN, this, 0, 0, VT, i, 0, (1.0 / S[i]), 0.0, 0, modeCount);
                        modeCount++;
                }
                t2 = MPI_Wtime();
                cout << "MOS SVD complete in " << t2 - t1 << " seconds" << endl;

                return 0;
        }

        else
        {
                cout << "min M not supported yet" << endl;
                throw(-1);
        }
}

int pMat::qr_run(int m, int n, int ia, int ja, std::vector<int> &ipiv)
{
        cout << "QR initializing" << endl;
        int IA = ia + 1;
        int JA = ja + 1;

        int JAnpm1 = JA + n - 1;
        int JAminMNm1 = JA + std::min(m, n) - 1;

        int ipiv_LOCc = std::max(1, numroc_(&JAnpm1, &nb, &(pG->mycol), &(desc[7]), &(pG->pcol)));
        int tau_LOCc = std::max(1, numroc_(&JAminMNm1, &nb, &(pG->mycol), &(desc[7]), &(pG->pcol)));
        cout << "ipiv_LOCc= " << ipiv_LOCc << endl;
        cout << "tau_LOCc= " << tau_LOCc << endl;
        vector<double> tau(1);
        vector<double> WORK(1);

        ipiv.resize(ipiv_LOCc);
        tau.resize(tau_LOCc);

        double t1, t2;
        int info = 0;
        int LWORK = -1;

        pdgeqpf(&m, &n, dataD.data(), &IA, &JA, desc, ipiv.data(), tau.data(), WORK.data(), &LWORK, &info);
        cout << "WORK=" << WORK[0] << " ,LWORK= " << LWORK << ", info= " << info << endl;
        LWORK = WORK[0];
        WORK.resize(LWORK);
        std::cout << "WORK allocated starting QR" << std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        t1 = MPI_Wtime();
        pdgeqpf(&m, &n, dataD.data(), &IA, &JA, desc, ipiv.data(), tau.data(), WORK.data(), &LWORK, &info);
        MPI_Barrier(MPI_COMM_WORLD);
        t2 = MPI_Wtime();
        cout << "QR complete in " << t2 - t1 << " seconds info =" << info << endl;
        MPI_Barrier(MPI_COMM_WORLD);
        tau.resize(0);
        WORK.resize(0);

        MPI_Barrier(MPI_COMM_WORLD);
        int ONE = 1;
        MPI_File fH;
        std::string pivot_name = "P.bin";
        MPI_File_open(MPI_COMM_WORLD, pivot_name.c_str(), MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &fH);
        if (printRank)
        {
                MPI_File_write(fH, &ONE, 1, MPI_INT, MPI_STATUS_IGNORE);
                MPI_File_write(fH, &N, 1, MPI_INT, MPI_STATUS_IGNORE);
                cout << "Write Start" << endl;
                cout << "M=" << ONE << "mb=" << mb << "N=" << N << "nb=" << nb << endl;
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
        if (printRank)
                cout << "MPI Allocation " << mpiEls << " , pblacs Allocation " << myRC[1] << endl;
        t1 = MPI_Wtime();
        MPI_File_set_view(fH, disp, MPI_INT, darray, "native", MPI_INFO_NULL);
        MPI_File_write_all(fH, ipiv.data(), mpiEls, MPI_INT, MPI_STATUS_IGNORE);
        MPI_File_close(&fH);
        t2 = MPI_Wtime();
        if (printRank)
                cout << "Write time is " << t2 - t1 << endl;

        return 1;
}

int pMat::transpose(pMat *A)
{
        transpose(A, this->M, this->N, 0, 0);
}

int pMat::transpose(pMat *A, int m, int n, int ia, int ja)
{
        int IA = ia + 1;
        int JA = ja + 1;
        cout << "Copying transpose" << endl;

        if ((m != M) && (n != A->M))
        {
                cout << "transpose dimension mismatch" << endl;
                cout << m << " " << M << " " << n << " " << A->M << endl;
                return -1;
        }
        int i_one = 1;
        double ONE = 1.0;
        double ZERO = 0.0;
        pdtran(&m, &n, &ONE, A->dataD.data(), &IA, &JA, A->desc, &ZERO, dataD.data(), &i_one, &i_one, desc);
}
int pMat::changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb, bool stdout)
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
                if(stdout)
                        cout << "Double changed pGrid" << endl;
                pdgemr2d(&m, &n, A->dataD.data(), &IA, &JA, A->desc, dataD.data(), &IB, &JB, desc, &(pG->icntxt));
        }
        if (type == 1)
        {
                if(stdout)
                        cout << "Complex changed pGrid" << endl;
                pzgemr2d(&m, &n, A->dataC.data(), &IA, &JA, A->desc, dataC.data(), &IB, &JB, desc, &(pG->icntxt));
        }
}
int pMat::changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb)
{
        changeContext(A, m, n, ia, ja, ib, jb, true);
}
int pMat::changeContext(pMat *A)
{
        changeContext(A, M, N, 0, 0, 0, 0);
}

int pMat::dMax(int dim, int rc, double &val)
{
        int index = 0;

        if (dim == 0)
        {
                int IA = 1, JA = rc + 1, i_one = 1;
                pdamax(&M, &val, &index, dataD.data(), &IA, &JA, desc, &i_one);
        }
        if (dim == 1)
        {
                int IA = rc + 1, JA = 1, i_one = 1;
                pdamax(&N, &val, &index, dataD.data(), &IA, &JA, desc, &M);
        }
}
int pMat::dSum(int dim, int rc, double &val)
{

        if (dim == 0)
        {
                int IA = 1, JA = rc + 1, i_one = 1;
                pdasum(&M, &val, dataD.data(), &IA, &JA, desc, &i_one);
        }
        if (dim == 1)
        {
                int IA = rc + 1, JA = 1, i_one = 1;
                pdasum(&N, &val, dataD.data(), &IA, &JA, desc, &M);
        }
}

void pMat::pinv(pMat *A)
{
        pMat *UU, *VV;
        UU = new pMat(A->M, std::min(A->M, A->N), A->pG);
        VV = new pMat(std::min(A->M, A->N), A->N, A->pG);
        vector<double> SS(std::min(A->M, A->N), 0.0);

        A->svd_run(A->M, A->N, 0, 0, UU, VV, SS);
        UU->write_bin("UU.bin");
        VV->write_bin("VV.bin");
        cout << "summing outer products" << endl;
        cout << "tol check " << SS[0] << " " << std::numeric_limits<double>::epsilon() * std::max(A->M, A->N) * SS[0] << "\r";
        this->matrix_Product('T', 'T', VV->N, UU->M, 1, VV, 0, 0, UU, 0, 0, 1.0 / SS[0], 0.0, 0, 0);
        for (int i = 1; i < SS.size(); i++)
        {
                //check tolerance
                cout << "tol check " << SS[i] << " " << std::numeric_limits<double>::epsilon() * std::max(A->M, A->N) * SS[0] << "\r";
                if (SS[i] > std::numeric_limits<double>::epsilon() * std::max(A->M, A->N) * SS[0])
                        this->matrix_Product('T', 'T', VV->N, UU->M, 1, VV, i, 0, UU, 0, i, 1.0 / SS[i], 1.0, 0, 0);
        }
        delete UU;
        delete VV;
}
int pMat::outerProductSum(pMat *U, char UT, pMat *VT, char VTT, std::vector<double> &S, int inv)
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
ostream &operator<<(std::ostream &os, const pMat &p)
{
        std::cout << "Descriptor type: " << p.desc[0] << std::endl;
        std::cout << "BLACS context: " << p.desc[1] << std::endl;
        std::cout << "Global Rows: " << p.desc[2] << std::endl;
        std::cout << "Global Cols: " << p.desc[3] << std::endl;
        std::cout << "Row Blocking factor: " << p.desc[4] << std::endl;
        std::cout << "Column Blocking factor: " << p.desc[5] << std::endl;
        std::cout << "Process row where first row is: " << p.desc[6] << std::endl;
        std::cout << "Process Col where first col is: " << p.desc[7] << std::endl;
        std::cout << "Leading Dimension: " << p.desc[8] << std::endl;
        std::cout << "Memory usage(data only) MB = " << p.MBs << std::endl;
        return os;
}
double pMat::getElement(int I, int J)
{
        double item = 0.0;
        int l, m;
        int x, y;

        l = I / (pG->prow * mb);
        m = J / (pG->pcol * nb);

        x = I % mb;
        y = J % nb;
        double temp = 0;
        if ((pG->myrow == (I / mb) % pG->prow) && (pG->mycol == (J / nb) % pG->pcol))
        {
                assert(((m * nb + y) * myRC[0] + l * mb + x) < nelements);
                temp = dataD[(m * nb + y) * myRC[0] + l * mb + x];
        }
        MPI_Request request;
        MPI_Iallreduce(MPI_IN_PLACE, &temp, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, &request);
        MPI_Wait(&request, MPI_STATUS_IGNORE);
        return temp;
}

double pMat::getLocalElement(int I, int J)
{
        double item = 0.0;
        int l, m;
        int x, y;

        l = I / (pG->prow * mb);
        m = J / (pG->pcol * nb);

        x = I % mb;
        y = J % nb;
        double temp = 0;
        if ((pG->myrow == (I / mb) % pG->prow) && (pG->mycol == (J / nb) % pG->pcol))
        {
                assert(((m * nb + y) * myRC[0] + l * mb + x) < nelements);
                temp = dataD[(m * nb + y) * myRC[0] + l * mb + x];
        }
        return temp;
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

int pMat::commCreate(MPI_Comm &col_comm, int dim)
{
        if (dim == 0)
        {
                cout << pG->mycol << " " << pG->rank << endl;
                MPI_Comm_split(MPI_COMM_WORLD, pG->mycol, pG->rank, &col_comm);
        }
        else if (dim == 1)
        {
                MPI_Comm_split(MPI_COMM_WORLD, pG->myrow, pG->rank, &col_comm);
        }
        else
                throw(-1);

        cout << col_comm << endl;
}
