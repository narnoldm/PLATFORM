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

    paramMap inputFile("QR_pre.inp", rank);

    double t2, t1;

    int debug_proc = 0;
    inputFile.getParamInt("stdout_proc", debug_proc);
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }

    string input;
    bool FOMInput = false;
    inputFile.getParamBool("BasisfromFOMData", FOMInput);
    if (FOMInput)
        inputFile.getParamString("inputStringFOM", input);

    bool BasisInput = false;
    inputFile.getParamBool("BasisfromDisk", BasisInput);
    if (BasisInput)
        inputFile.getParamString("inputStringBasis", input);

    if (BasisInput == FOMInput)
    {
        cout << "can't have both inputs or no inputs" << endl;
        return 0;
    }

    int numModes;
    inputFile.getParamInt("numModes", numModes);
    double pSampling;
    inputFile.getParamDouble("pSampling", pSampling);
    int numModesUsol = numModes;

    //
    string deimFolder, deimPrefix, dfd_file;
    inputFile.getParamString("deimModesFolder", deimFolder);
    inputFile.getParamString("deimModesPrefix", deimPrefix);
    inputFile.getParamString("dfd_itype_file", dfd_file);

    bool loadRHS=false;
    string input2;
    inputFile.getParamBool("RHS load",loadRHS);
    if(loadRHS)
    {
        inputFile.getParamString("inputStringRHS",input2);
    }


    cout << "input string is: " << input << endl;
    vector<string> token;
    tokenparse(input, "|", token);
    meta *dataset1,*dataset2;
    string firstFile;
    if (FOMInput)
    {
        dataset1 = new tecIO(token);
        firstFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;
        dynamic_cast<tecIO *>(dataset1)->activateReorder(firstFile.c_str());
    }

    if (BasisInput)
    {
        dataset1 = new meta(token);
        firstFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;
    }

    PGrid *evenG;
    evenG = new PGrid(rank, size, 0);

    pMat *A, *U, *VT;
    
    if (FOMInput)
    {
        A = new pMat(dataset1->nPoints, dataset1->nSets, evenG, 0, 0, 0.0);
        t1 = MPI_Wtime();
        dataset1->batchRead(A);
        t2 = MPI_Wtime();
        cout << "load took " << t2 - t1 << " seconds" << endl;

        t1 = MPI_Wtime();

        tecIO *datasetTec = dynamic_cast<tecIO *>(dataset1);
        datasetTec->calcAvg(A);
        datasetTec->subAvg(A);
        datasetTec->calcNorm(A);
        datasetTec->normalize(A);

        t2 = MPI_Wtime();
        cout << "preprocessing took " << t2 - t1 << " seconds" << endl;
        //A->write_bin("A.bin");
        int M = dataset1->nPoints, N = dataset1->nSets;

        vector<double> S;

        U = new pMat(M, std::min(M, N), evenG);
        VT = new pMat(std::min(M, N), N, evenG);
        S.resize(std::min(M, N));

        A->svd_run(M, N, 0, 0, U, VT, S);
        delete A;
        delete VT;
    }
    int PointsNeeded = 0;
    if (FOMInput)
        PointsNeeded = dynamic_cast<tecIO *>(dataset1)->nCells * pSampling;

    int nCells = 0, nVars = 0;
    if (BasisInput)
    {
        inputFile.getParamInt("nCells", nCells);
        inputFile.getParamInt("nVars", nVars);
        PointsNeeded = nCells * pSampling;
    }
    vector<int> P;
    pMat *UsT, *Uread;
    if (FOMInput)
    {
        Uread = new pMat(U->M, numModes, evenG);
        UsT = new pMat(numModes, U->M, evenG);
        Uread->changeContext(U, U->M, numModes, 0, 0, 0, 0);
        UsT->transpose(U, UsT->M, UsT->N, 0, 0);
        delete U;
        Uread;
    }
    if (BasisInput)
    {
        Uread = new pMat(dataset1->nPoints, dataset1->nSets, evenG);
        dataset1->batchRead(Uread);
        UsT = new pMat(numModes, Uread->M, evenG);
        UsT->transpose(Uread);
        assert(UsT->N == nCells * nVars);
    }

    UsT->qr_run(UsT->M, UsT->N, 0, 0, P);

    vector<int> gP;
    vector<int> itype;
    set<int> samplingPoints;

    t1 = MPI_Wtime();
    if (rank == 0)
    {
        readMat("P.bin", gP);
        gP.resize(numModes);
        for (int i = 0; i < gP.size(); i++)
        {
            gP[i]--; //switch to 0 indexing
            //switch to physical points
            cout << i << " " << gP[i] << " ";
            if (FOMInput)
                gP[i] = gP[i] % dynamic_cast<tecIO *>(dataset1)->nCells;
            if (BasisInput)
                gP[i] = gP[i] % nCells;
            auto check = samplingPoints.emplace(gP[i]);
            if (!check.second)
            {
                cout << "repeated element" << endl;
            }
        }
        cout << "goal is " << PointsNeeded << "points" << endl;
        cout << "points after qr: " << samplingPoints.size() << " of " << PointsNeeded << endl;
        /*for (std::set<int>::iterator it = samplingPoints.begin(); it != samplingPoints.end(); ++it)
        {
            cout << *it << endl;
        }*/
        //add boundary cells
        readMat(dfd_file, itype);
        for (vector<int>::iterator it = itype.begin(); it != itype.end(); ++it)
        {
            if (*it != 0)
            {
                auto check = samplingPoints.emplace(it - itype.begin());
                if (!check.second)
                {
                    cout << "repeated element " << it - itype.begin() << "\r";
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
            vector<int> rPoints;
            if (FOMInput)
                rPoints.resize(dynamic_cast<tecIO *>(dataset1)->nCells, 0);
            if (BasisInput)
                rPoints.resize(nCells);
            for (int i = 0; i < rPoints.size(); i++)
                rPoints[i] = i;
            srand(1);
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
        vector<double> gPD;
        if (FOMInput)
            gPD.resize(dynamic_cast<tecIO *>(dataset1)->nCells, 0.0);
        if (BasisInput)
            gPD.resize(nCells, 0.0);
        std::copy(samplingPoints.begin(), samplingPoints.end(), gP.begin());

        for_each(gP.begin(), gP.end(), [](int &tt) { tt += 1; });
        printASCIIVecP0("samplingPoints.txt", gP, gP.size());
        writeMat("Pall.bin", gP.size(), 1, gP);
        for_each(gP.begin(), gP.end(), [](int &tt) { tt -= 1; });

        for (int i = 0; i < samplingPoints.size(); i++)
            gPD[gP[i]] = 1.0;

        vector<string> Pname;
        std::string tempname = "sampling";
        Pname.push_back(tempname);
        if (FOMInput)
            dynamic_cast<tecIO *>(dataset1)->writeSingleFile("sampling.szplt", Pname, gPD.data(), firstFile);

        if (BasisInput)
        {
            tecIO *tempPoint = new tecIO();
        }
    }

    t2 = MPI_Wtime();
    cout << "figuring out samples took " << t2 - t1 << " seconds" << endl;

    if (rank != 0)
    {
        gP.resize(PointsNeeded, 0);
    }
    cout << "broadcasting points" << endl;
    MPI_Bcast(gP.data(), gP.size(), MPI_INT, 0, MPI_COMM_WORLD);

    cout << "calculating DIEM interpolant" << endl;
    pMat *Usamp;

    if (FOMInput)
        Usamp = new pMat(gP.size() * dynamic_cast<tecIO *>(dataset1)->numVars, numModes, evenG);

    if (BasisInput)
        Usamp = new pMat(gP.size() * nVars, numModes, evenG);

    t1 = MPI_Wtime();
    for (int i = 0; i < gP.size(); i++)
    {
        cout << (double)i / gP.size() * 100 << "percent points extracted \r";
        if (FOMInput)
        {
            for (int j = 0; j < dynamic_cast<tecIO *>(dataset1)->numVars; j++)
                Usamp->changeContext(Uread, 1, numModes, gP[i] + j * dynamic_cast<tecIO *>(dataset1)->nCells, 0, i + j * gP.size(), 0, false);
        }
        if (BasisInput)
        {
            for (int j = 0; j < nVars; j++)
                Usamp->changeContext(Uread, 1, numModes, gP[i] + j * nCells, 0, i + j * gP.size(), 0, false);
        }
    }

    t2 = MPI_Wtime();
    cout << "extraction took " << t2 - t1 << " seconds" << endl;

    Usamp->write_bin("Usamp.bin");

    t1 = MPI_Wtime();
    pMat *pinvUsamp = new pMat(Usamp->N, Usamp->M, evenG);
    pinvUsamp->pinv(Usamp);

    pinvUsamp->write_bin("pinvUsamp.bin");

    //Assuming Usol=U and numModes= numModesSol This will reduce to I, but doing math for when we have RHS
    pMat *Usol = Uread;


    pMat *UsolU = new pMat(numModesUsol, numModes, evenG);

    UsolU->matrix_Product('T', 'N', numModesUsol, numModes, Uread->M, Usol, 0, 0, Uread, 0, 0, 1.0, 0.0, 0, 0);

    pMat *deimInterp = new pMat(UsolU->N, pinvUsamp->N, evenG);

    deimInterp->matrix_Product('N', 'N', deimInterp->M, deimInterp->N, pinvUsamp->M, UsolU, 0, 0, pinvUsamp, 0, 0, 1.0, 0.0, 0, 0);
    t2 = MPI_Wtime();
    cout << "DEIM interpolant calculation took " << t2 - t1 << " seconds" << endl;

    //deimInterp->write_bin("deimInterp.bin");

    //delete all the extras
    delete UsolU;
    delete pinvUsamp;
    delete Usamp;
    delete UsT;
    delete Uread;

    pMat *deimInterp_T = new pMat(deimInterp->N, deimInterp->M, deimInterp->pG);
    deimInterp_T->transpose(deimInterp);
    delete deimInterp;

    Uread = new pMat(dataset1->nPoints, numModes, evenG);
    // 0 insertion for output

    t1 = MPI_Wtime();
    for (int i = 0; i < gP.size(); i++)
    {
        cout << (double)i / gP.size() * 100 << "percent points emplaced \r";
        if (FOMInput)
        {
            for (int j = 0; j < dynamic_cast<tecIO *>(dataset1)->numVars; j++)
                Uread->changeContext(deimInterp_T, 1, numModes, i + j * gP.size(), 0, gP[i] + j * dynamic_cast<tecIO *>(dataset1)->nCells, 0, false);
        }
        if (BasisInput)
        {
            for (int j = 0; j < nVars; j++)
                Uread->changeContext(deimInterp_T, 1, numModes, i + j * gP.size(), 0, gP[i] + j * nCells, 0, false);
        }
    }
    t2 = MPI_Wtime();
    cout << "DEIM emplacement took " << t2 - t1 << " seconds" << endl;

    if (FOMInput)
    {
        tecIO *Uout = new tecIO();
        Uout->snap0 = 1;
        Uout->snapF = Uread->N;
        Uout->snapSkip = 1;
        Uout->nSets = Uread->N;
        //Uout->prefix = "deimInterp";
        //Uout->suffix = ".szplt";
        Uout->isInit = true;
        Uout->meshFile = dataset1->prefix + std::to_string(dataset1->snap0) + dataset1->suffix;
        Uout->fixedMesh = true;
        Uout->getDimNodes();
        Uout->varName = dynamic_cast<tecIO *>(dataset1)->varName;
        Uout->varIndex = dynamic_cast<tecIO *>(dataset1)->varIndex;
        Uout->numVars = Uout->varName.size();
        Uout->nPoints = Uout->nCells * Uout->numVars;

        Uout->activateReorder(firstFile.c_str());
        Uout->activateGEMSbin(firstFile.c_str());
        t1 = MPI_Wtime();
        Uout->batchWrite(Uread, deimFolder, deimPrefix);
        t2 = MPI_Wtime();
        cout << "Output took " << t2 - t1 << " seconds" << endl;
    }

    if (BasisInput)
    {
        meta *Uout = new meta();
        Uout->snap0 = 1;
        Uout->snapF = Uread->N;
        Uout->snapSkip = 1;
        Uout->nSets = Uread->N;
        Uout->prefix = "deimInterp";
        Uout->suffix = ".bin";
        Uout->isInit = true;

        Uout->nPoints = nCells * nVars;
        assert(Uout->nPoints == nCells * nVars);
        t1 = MPI_Wtime();
        Uout->batchWrite(Uread, deimFolder, deimPrefix);
        t2 = MPI_Wtime();
        cout << "Output took " << t2 - t1 << " seconds" << endl;
    }

    cout.rdbuf(strm_buffer);
    MPI_Finalize();

    return 0;
}
