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
        cout << "local processor row"<<myrow<<" and column "<<mycol<<endl;
        cout << "PBLACS is row major by default and everything is else is Column major (because reasons)"<<endl;
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
                mb = nb;
                        cout << "mb/nb = " << nb << endl;
                MPI_Barrier(MPI_COMM_WORLD);
        }
        else if (block == 1) //load blocks
        {
                nb = 128;
                mb = M;
                if (nb > (N / pG->getDim(0)))
                        nb = std::max(1,N / pG->getDim(0));
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
        else if (block ==3) //synched block
        {
                nb=N/(pG->prow); 
                mb=M/(pG->pcol);
        }
        myRC[0] =numroc_(&M, &mb, &(pG->myrow), &i_zero, &(pG->prow));
        //cout<<myRC[0]<<" "<<M<<" "<<mb<<" "<<pG->myrow<<" "<<pG->prow<<endl;
        
        myRC[1] = numroc_(&N, &nb, &(pG->mycol), &i_zero, &(pG->pcol));
        //cout<<myRC[1]<<" "<<N<<" "<<nb<<" "<<pG->mycol<<" "<<pG->pcol<<endl;

        //if (myRC[0] < 1)
        //        myRC[0] = 1;
        //if (myRC[1] < 1)
        //        myRC[1] = 1;
        nelements = (long long)myRC[0] * (long long)myRC[1];
        for(int i=0;i<2;i++)
                myRC[i]=std::max(1,myRC[i]);

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
        cout<<"nelements "<<nelements<<endl;
        cout<<myRC[0]<<" "<<myRC[1]<<endl;

        

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
                                        cout <<i<<" "<< dataD[i] << endl;
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
                cout << "Write Start" << endl;
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
        cout<<"M = "<<M<<" N = "<<N<<endl;
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
        int info=0;
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

        if(info<0)
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
int pMat::mos_run(int M, int N, int ia, int ja, pMat *&U, pMat *&VT, vector<double> &S)
{
        int IA = ia + 1;
        int JA = ja + 1;
        int i_one = 1;
        double t2, t1;
        cout<<"starting MOS"<<endl;
        int minMN=std::min(M,N);
        pMat *corMat = new pMat(minMN,minMN,pG,0,0,0.0);
        pMat * corMatp0 = new pMat(minMN,minMN,pG,0,2,0.0);
        pMat * VTp0 = new pMat(minMN,minMN,pG,0,2,0.0);
        if(minMN==N)
        {
                corMat->matrix_Product('T','N',minMN,minMN,M,this,0,0,this,0,0,1.0,0.0,0,0);
                cout<<"cor Mat created"<<endl;
                corMatp0->changeContext(corMat);
                cout<<"start eigensolve"<<endl;
                 t1 = MPI_Wtime();
                if(pG->rank==0)
                {
                Eigen::MatrixXd corMatEig= Eigen::Map<Eigen::MatrixXd>(corMatp0->dataD.data(),minMN,minMN);
                corMatEig = corMatEig.selfadjointView<Eigen::Upper>();
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(corMatEig);
                
                Eigen::MatrixXd VTEig(minMN,minMN);
                Eigen::VectorXd SEig(minMN);

                SEig=es.eigenvalues().reverse();
                VTEig=es.eigenvectors().rowwise().reverse();
                

                Eigen::Map<Eigen::MatrixXd>(VTp0->dataD.data(),VTEig.rows(),VTEig.cols())=VTEig;
                Eigen::Map<Eigen::VectorXd>(S.data(),minMN)=SEig;

                for(int i=0;i<S.size();i++)
                        S[i]=std::sqrt(S[i]);


                }
                cout<<"finish eigensolve"<<endl;
                MPI_Bcast(S.data(), S.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
                pMat *V = new pMat(VT->M,VT->N,VT->pG,0,VT->block,0.0);
                V->changeContext(VTp0);
                VT->transpose(V);
                delete V;

                delete VTp0;
                delete corMatp0;

                for(int i=0;i<minMN;i++)
                {
                        cout<<i<<endl;
                        U->matrix_Product('N','T',M,1,minMN,this,0,0,VT,i,0,(1.0/S[i]),0.0,0,i);
                }
                t2 = MPI_Wtime();
                cout << "MOS SVD complete in " << t2 - t1 << " seconds" << endl;

        }
        


        else
        {
                cout<<"min M not supported yet"<<endl;
                throw(-1);
        }
}

int pMat::transpose(pMat *A)
{
        transpose(A,this->M,this->N,0,0);
}

int pMat::transpose(pMat *A, int m, int n, int ia, int ja)
{
        int IA = ia + 1;
        int JA = ja + 1;
        if (printRank)
                cout << "Copying transpose" << endl;

        if ((m != M) && (n != A->M))
                if (printRank)
                {
                        cout << "transpose dimension mismatch" << endl;
                        return -1;
                }
        int i_one = 1;
        double ONE = 1.0;
        double ZERO = 0.0;
        pdtran(&m, &n, &ONE, A->dataD.data(), &IA, &JA, A->desc, &ZERO, dataD.data(), &i_one, &i_one, desc);
}
int pMat::changeContext(pMat *A, int m, int n, int ia, int ja, int ib, int jb)
{
        int IA = ia + 1;
        int JA = ja + 1;
        int IB = ib + 1;
        int JB = jb + 1;
                cout << "Copying Matrix" << endl
                     << "m = " << m << " , "
                     << "n = " << n << endl;
        int i_one = 1;
        if (type == 0)
        {
                        cout << "Double changed pGrid" << endl;
                pdgemr2d(&m, &n, A->dataD.data(), &IA, &JA, A->desc, dataD.data(), &IB, &JB, desc, &(pG->icntxt));
        }
        if (type == 1)
        {
                        cout << "Complex changed pGrid" << endl;
                pzgemr2d(&m, &n, A->dataC.data(), &IA, &JA, A->desc, dataC.data(), &IB, &JB, desc, &(pG->icntxt));
        }
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


int pMat::outerProductSum(pMat *U,char UT,pMat *VT,char VTT,std::vector<double> &S,int inv)
{
        if(inv==1)
        {
                this->matrix_Product(UT,VTT,U->M,VT->N,1,U,0,0,VT,0,0,1/S[0],0.0,0,0);
                for(int i=1;i<S.size();i++)
                {
                        this->matrix_Product(UT,VTT,U->M,VT->N,1,U,0,i,VT,i,0,1/S[i],1.0,0,0);
                }
        }
        else
        {
                this->matrix_Product(UT,VTT,U->M,VT->N,1,U,0,0,VT,0,0,S[0],0.0,0,0);
                for(int i=1;i<S.size();i++)
                {
                        this->matrix_Product(UT,VTT,U->M,VT->N,1,U,0,i,VT,i,0,S[i],1.0,0,0);       
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
        int x,y;

        l= I/(pG->prow*mb);
        m=J/(pG->pcol*nb);

        x= I%mb;
        y= J%nb;
        double temp=0;
        if((pG->myrow== (I/mb)%pG->prow) && (pG->mycol== (J/nb)%pG->pcol))
        {
                assert(((m*nb+y)*myRC[0]+ l*mb +x)<nelements);
                temp=dataD[(m*nb+y)*myRC[0]+ l*mb +x];
        }
        MPI_Request request;
        MPI_Iallreduce(MPI_IN_PLACE,&temp,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,&request);
        MPI_Wait(&request,MPI_STATUS_IGNORE);
        return temp;
}

double pMat::getLocalElement(int I, int J)
{
        double item = 0.0;
        int l, m;
        int x,y;

        l= I/(pG->prow*mb);
        m=J/(pG->pcol*nb);

        x= I%mb;
        y= J%nb;
        double temp=0;
        if((pG->myrow== (I/mb)%pG->prow) && (pG->mycol== (J/nb)%pG->pcol))
        {
                assert(((m*nb+y)*myRC[0]+ l*mb +x)<nelements);
                temp=dataD[(m*nb+y)*myRC[0]+ l*mb +x];
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
                                double epsilon=.01;
                                if(p1.dataD[i]>0)
                                {
                                if ( (p1.dataD[i]+p1.dataD[i]*epsilon) < p2.dataD[i] || (p1.dataD[i]-p1.dataD[i]*epsilon) > p2.dataD[i])
                                {
                                        cout << "element " << i << " does not match" << endl;
                                        cout<<(p1.dataD[i]-p1.dataD[i]*epsilon)<<"< "<< p2.dataD[i]<<" < "<<(p1.dataD[i]+p1.dataD[i]*epsilon)<<endl;
                                        cout<<p1.dataD[i]<<"!="<<p2.dataD[i]<<endl;
                                        return false;
                                }
                                }
                                else
                                {
                                        if ( (p1.dataD[i]+p1.dataD[i]*epsilon) > p2.dataD[i] || (p1.dataD[i]-p1.dataD[i]*epsilon) < p2.dataD[i])
                                {
                                        cout << "element " << i << " does not match" << endl;
                                        cout<<(p1.dataD[i]-p1.dataD[i]*epsilon)<<"< "<< p2.dataD[i]<<" < "<<(p1.dataD[i]+p1.dataD[i]*epsilon)<<endl;
                                        cout<<p1.dataD[i]<<"!="<<p2.dataD[i]<<endl;
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

int pMat::commCreate(MPI_Comm &col_comm,int dim)
{
        if(dim==0)
        {
                cout<<pG->mycol<<" "<<pG->rank<<endl;
                MPI_Comm_split(MPI_COMM_WORLD,pG->mycol,pG->rank,&col_comm);
        }
        else if(dim==1)
        {
                        MPI_Comm_split(MPI_COMM_WORLD,pG->myrow,pG->rank,&col_comm);
        }
        else
                throw(-1);


        cout<<col_comm<<endl;
}
