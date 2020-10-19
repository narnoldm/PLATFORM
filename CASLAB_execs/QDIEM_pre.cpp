#include "metadata.hpp"
#include "param.hpp"

#include <set>

using namespace ::std;
int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();

    paramMap inputFile("QR_pre.inp",rank);


    int debug_proc;
    inputFile.getParamInt("stdout_proc",debug_proc);
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }


    string input;
    inputFile.getParamString("inputString",input);


    cout << "input string is: " << input << endl;
    vector<string> token;
    tokenparse(input, "|", token);

    tecIO *dataset1 = new tecIO(token);
    string firstFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;
    dataset1->activateReorder(firstFile.c_str());

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);

    pMat *A;

    A = new pMat(dataset1->nPoints, dataset1->nSets, evenG, 0, 0, 0.0);
    dataset1->batchRead(A);
    dataset1->calcAvg(A);
    dataset1->subAvg(A);
    dataset1->calcNorm(A);
    dataset1->normalize(A);

    A->write_bin("A.bin");
    int M = dataset1->nPoints, N = dataset1->nSets;

    int numModes = std::min(M, N) * .8;
    double pSampling = .2;
    int numModesUsol = numModes;

    int PointsNeeded = dataset1->nCells * pSampling;

    pMat *U, *VT, *UsT;
    vector<double> S;

    U = new pMat(M, std::min(M, N), evenG);
    VT = new pMat(std::min(M, N), N, evenG);
    S.resize(std::min(M, N));

    A->svd_run(M, N, 0, 0, U, VT, S);
    for (int i = 0; i < numModes; i++)
        cout << S[i] << endl;
    vector<int> P;

    UsT = new pMat(numModes, U->M, evenG);
    UsT->transpose(U, UsT->M, UsT->N, 0, 0);

    UsT->qr_run(UsT->M, UsT->N, 0, 0, P);

    vector<int> gP;
    vector<int> itype;
    set<int> samplingPoints;
    if (rank == 0)
    {
        readMat("P.bin", gP);
        gP.resize(numModes);
        for (int i = 0; i < gP.size(); i++)
        {
            gP[i]--; //switch to 0 indexing
            //switch to physical points
            cout << i << " " << gP[i] << " ";
            gP[i] = gP[i] % dataset1->nCells;
            cout << gP[i] << endl
                 << endl;
            auto check = samplingPoints.emplace(gP[i]);
            if (!check.second)
            {
                cout << "repeated element" << endl;
            }
        }
        cout << "goal is " << PointsNeeded << "points" << endl;
        cout << "points after qr: " << samplingPoints.size() << " of " << PointsNeeded << endl;
        for (std::set<int>::iterator it = samplingPoints.begin(); it != samplingPoints.end(); ++it)
        {
            cout << *it << endl;
        }
        //add boundary cells
        readMat("dfd_itype.bin", itype);
        for (vector<int>::iterator it = itype.begin(); it != itype.end(); ++it)
        {
            if (*it != 0)
            {
                auto check = samplingPoints.emplace(it - itype.begin());
                if (!check.second)
                {
                    cout << "repeated element " << it - itype.begin() <<"\r";
                }
            }
        }
        cout << "points after boundaries: " << samplingPoints.size() << " of " << PointsNeeded << endl;
        assert(PointsNeeded > samplingPoints.size());
        cout << "Need " << PointsNeeded - samplingPoints.size() << " more points" << endl;
        int r = 0;
        while (PointsNeeded - samplingPoints.size() > 0)
        {
            //Random Points
            vector<int> rPoints(dataset1->nCells, 0);
            for (int i = 0; i < rPoints.size(); i++)
                rPoints[i] = i;

            random_shuffle(rPoints.begin(), rPoints.end());
            rPoints.resize(PointsNeeded - samplingPoints.size());
            for (vector<int>::iterator it = rPoints.begin(); it != rPoints.end(); ++it)
            {
                if (*it != 0)
                {
                    auto check = samplingPoints.emplace(*it);
                    if (!check.second)
                    {
                        cout << "repeated element " << *it << "\r";
                    }
                }
            }
            r++;
            cout << "after sweep " << r << " need " << PointsNeeded - samplingPoints.size() << " more Points" << endl;
            if (r >= 20)
            {
                cout << "still don't have the points" << endl;
                throw(-1);
            }
        }
        cout << "All points found " << endl;
        gP.resize(samplingPoints.size(), 0);
        vector<double> gPD(dataset1->nCells, 0.0);
        std::copy(samplingPoints.begin(), samplingPoints.end(), gP.begin());
        printASCIIVecP0("samplingPoints.txt", gP, gP.size());
        writeMat("Pall.bin", gP.size(), 1, gP);
        for (int i = 0; i < samplingPoints.size(); i++)
            gPD[gP[i]] = 1.0;

        vector<string> Pname;
        std::string tempname = "sampling";
        Pname.push_back(tempname);
        dataset1->writeSingleFile("sampling.szplt", Pname, gPD.data(), firstFile);
    }
    if (rank != 0)
    {
        gP.resize(PointsNeeded, 0);
    }
    cout << "broadcasting points" << endl;
    MPI_Bcast(gP.data(), gP.size(), MPI_INT, 0, MPI_COMM_WORLD);

    cout << "calculating DIEM interpolant" << endl;

    pMat *Usamp = new pMat(gP.size() * dataset1->numVars, numModes, evenG);
    for (int i = 0; i < gP.size(); i++)
    {
        cout << i << "\r";
        for (int j = 0; j < dataset1->numVars; j++)
            Usamp->changeContext(U, 1, numModes, gP[i] + j * dataset1->nCells, 0, i + j * gP.size(), 0,false);
    }

    Usamp->write_bin("Usamp.bin");

    pMat *pinvUsamp = new pMat(Usamp->N, Usamp->M, evenG);
    pinvUsamp->pinv(Usamp);

    pinvUsamp->write_bin("pinvUsamp.bin");

    //Assuming Usol=U and numModes= numModesSol This will reduce to I, but doing math for when we have RHS
    pMat *Usol = U;
    pMat *UsolU = new pMat(numModesUsol, numModes, evenG);

    UsolU->matrix_Product('T', 'N', numModesUsol, numModes, U->M, Usol, 0, 0, U, 0, 0, 1.0, 0.0, 0, 0);

    pMat *deimInterp = new pMat(UsolU->N, pinvUsamp->N, evenG);

    deimInterp->matrix_Product('N', 'N', deimInterp->M, deimInterp->N, pinvUsamp->M, UsolU, 0, 0, pinvUsamp, 0, 0, 1.0, 0.0, 0, 0);

    deimInterp->write_bin("deimInterp.bin");

    //delete all the extras
    delete UsolU;
    delete pinvUsamp;
    delete Usamp;
    delete UsT;

    pMat *deimInterp_T = new pMat(deimInterp->N, deimInterp->M, deimInterp->pG);
    deimInterp_T->transpose(deimInterp);
    delete deimInterp;
    delete U;

    U = new pMat(M, numModes, evenG);
    // 0 insertion for output
    for (int i = 0; i < gP.size(); i++)
    {
        cout << i <<"\r";
        for (int j = 0; j < dataset1->numVars; j++)
            U->changeContext(deimInterp_T, 1, numModes, i + j * gP.size(), 0, gP[i] + j * dataset1->nCells, 0,false);
    }

    tecIO *Uout = new tecIO();
    Uout->snap0 = 1;
    Uout->snapF = U->N;
    Uout->snapSkip = 1;
    Uout->nSets = U->N;
    Uout->prefix = "deimInterp";
    Uout->suffix = ".szplt";
    Uout->isInit = true;
    Uout->meshFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;
    Uout->fixedMesh = true;
    Uout->getDimNodes();
    Uout->varName = dataset1->varName;
    Uout->varIndex = dataset1->varIndex;
    Uout->numVars = Uout->varName.size();
    Uout->nPoints = Uout->nCells * Uout->numVars;

    Uout->activateReorder(firstFile.c_str());
    Uout->activateGEMSbin(firstFile.c_str());
    Uout->batchWrite(U, "deimInterp", "Deim_mode_");

    /*string filename = "gemsma1.bin";

    if (A->check_bin_size(filename, M, N))
    {

        A = new pMat(M, N, evenG);
        A->read_bin(filename);

        pMat *U, *VT, *UsT;
        vector<double> S;

        U = new pMat(M, std::min(M, N), evenG);
        VT = new pMat(std::min(M, N), N, evenG);
        S.resize(std::min(M, N));

        A->svd_run(M, N, 0, 0, U, VT, S);
        for (int i = 0; i < numModes; i++)
            cout << S[i] << endl;
        vector<int> P;

        UsT = new pMat(numModes, U->M, evenG);
        UsT->transpose(U, UsT->M, UsT->N, 0, 0);

        UsT->qr_run(UsT->M, UsT->N, 0, 0, P);

        vector<int> gP;
        if (!rank)
        {
            readMat("P.bin", gP);
            gP.resize(numModes);
            for (int i = 0; i < gP.size(); i++)
            {
                gP[i]--; //switch to 0 indexing
                cout << i << " " << gP[i] << endl;
                //switch to physical points
            }
            //add boundary points

            //remove duplicate physical points
        }
    }*/

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}
