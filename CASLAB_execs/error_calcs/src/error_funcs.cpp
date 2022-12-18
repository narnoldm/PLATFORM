#include "error_funcs.hpp"

using namespace :: std;

void write_failed_abs_and_l2_error(tecIO* setData, string errDir, string errSuffix, bool mags, vector<string> magsToken)
{
    // do some pre-processing if calculating magnitudes
    int numVars;
    vector<int> magFactor_int;
    int nUniqueGroups;
    vector<vector<int>> groupRef;
    vector<string> varNames;
    if (mags)
    {

        // get sorted list of unique norm factors
        magFactor_int.resize((setData->numVars));
        for (int i = 0; i < setData->numVars; ++i)
        {
            magFactor_int[i] = stoi(magsToken[i]);
            if (i > 0)
            {
                if (magFactor_int[i] < magFactor_int[i-1])
                {
                    cout << "magsInput must be strictly non-decreasing" << endl;
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }
        }

        std::sort(magFactor_int.begin(), magFactor_int.end());                  // sort
        auto last = std::unique(magFactor_int.begin(), magFactor_int.end());    // get unique elements
        magFactor_int.erase(last, magFactor_int.end());

        nUniqueGroups = magFactor_int.size();

        groupRef.resize(nUniqueGroups);
        varNames.resize(nUniqueGroups);
        numVars = nUniqueGroups;

        // get indices within dataset for variables in each group
        for (int i = 0; i < setData->numVars; ++i)
        {
            for (int j = 0; j < nUniqueGroups; ++j)
            {
                if (stoi(magsToken[i]) == magFactor_int[j])
                {
                    groupRef[j].push_back(i);
                    varNames[j] += setData->varName[i];
                }
            }
        }

        for (int j = 0; j < nUniqueGroups; ++j)
        {
            if (groupRef[j].size() > 1)
            {
                varNames[j] += "Mag";
            }
        }

        if (numVars == setData->numVars)
        {
            cout << "Requested magnitude calculations, but no matching variables" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        numVars = setData->numVars;
        varNames = setData->varName;
    }

    int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        ofstream out;
        string outFile;
        outFile = errDir + "/abs_avg_sum_err" + errSuffix + ".dat";
		out.open(outFile, ios::trunc);
		for (int i = 0; i < numVars; ++i) {
            out << varNames[i] + ": nan" << endl;
        }
        out << "Average: nan" << endl;
        out.close();

        outFile = errDir + "/abs_rel_sum_err" + errSuffix + ".dat";
		out.open(outFile, ios::trunc);
		for (int i = 0; i < numVars; ++i) {
            out << varNames[i] + ": nan" << endl;
        }
        out << "Average: nan" << endl;
        out.close();

        outFile = errDir + "/l2_avg_sum_err" + errSuffix + ".dat";
		out.open(outFile, ios::trunc);
		for (int i = 0; i < numVars; ++i) {
            out << varNames[i] + ": nan" << endl;
        }
        out << "Average: nan" << endl;
        out.close();

        outFile = errDir + "/l2_rel_sum_err" + errSuffix + ".dat";
		out.open(outFile, ios::trunc);
		for (int i = 0; i < numVars; ++i) {
            out << varNames[i] + ": nan" << endl;
        }
        out << "Average: nan" << endl;
        out.close();

    }

}

void calc_integrated_error(pMat* dataMat, vector<string>& varNames, string outFile) {

    int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    pMat* onesCol = new pMat(dataMat->N, 1, dataMat->pG, 0, 0, 1.0, false);
    pMat* dataMatInt = new pMat(dataMat->M, 1, dataMat->pG, false);
    pMat* dataMatIntP0 = new pMat(dataMat->M, 1, dataMat->pG, 0, 2, 0.0, false);
    double dataMatIntAvg = 0.0;

    dataMatInt->matrix_Product('N', 'N', dataMat->M, 1, dataMat->N, dataMat, 0, 0, onesCol, 0, 0, 1.0 / dataMat->N, 0.0, 0, 0);
    dataMatIntP0->changeContext(dataMatInt, false);
    if (rank == 0) {
        ofstream out;
		out.open(outFile, ios::trunc);
		for (int i = 0; i < dataMatIntP0->dataD.size(); ++i) {
            out << varNames[i] + ": " << setprecision(numeric_limits<double>::digits10) << dataMatIntP0->dataD[i] << endl;
            dataMatIntAvg += dataMatIntP0->dataD[i];
        }
        dataMatIntAvg /= dataMat->M;
        out << "Average: " << setprecision(numeric_limits<double>::digits10) << dataMatIntAvg << endl;
        out.close();
    }

    destroyPMat(onesCol, false);
    destroyPMat(dataMatInt, false);
    destroyPMat(dataMatIntP0, false);

}

void calc_magnitude_data(pMat* dataMat, pMat* dataMatMags, tecIO* setData, vector<vector<int>> groupRef)
{
    // NOTE: this assumes that the data has been reordered already

    // square
    for (int j = 0; j < dataMat->dataD.size(); ++j)
    {
        dataMat->dataD[j] = dataMat->dataD[j] * dataMat->dataD[j];
    }


    // square and add for each group
    for (int j = 0; j < groupRef.size(); ++j)
    {
        for (int k = 0; k < groupRef[j].size(); ++k)
        {
            dataMatMags->matrix_Sum('N', setData->nCells, setData->nSets, dataMat, groupRef[j][k] * setData->nCells, 0, j * setData->nCells, 0, 1.0, 1.0);
        }
    }

    // square root
    for (int j = 0; j < dataMatMags->dataD.size(); ++j)
    {
        dataMatMags->dataD[j] = sqrt(dataMatMags->dataD[j]);
    }

}

void calc_abs_and_l2_error(pMat* dataTruth, pMat* dataComp, tecIO* setData, string errDir, string errSuffix, bool outErrField)
{
    vector<string> magsToken = {""};
    calc_abs_and_l2_error(dataTruth, dataComp, setData, errDir, errSuffix, outErrField, false, magsToken);
}

void calc_abs_and_l2_error(pMat* dataTruth, pMat* dataComp, tecIO* setData, string errDir, string errSuffix, bool outErrField,
        bool mags, vector<string> magsToken)
{

    // do some pre-processing if calculating magnitudes
    int numVars;
    vector<int> magFactor_int;
    int nUniqueGroups;
    vector<vector<int>> groupRef;
    vector<string> varNames;
    if (mags)
    {

        // get sorted list of unique norm factors
        magFactor_int.resize((setData->numVars));
        for (int i = 0; i < setData->numVars; ++i)
        {
            magFactor_int[i] = stoi(magsToken[i]);
            if (i > 0)
            {
                if (magFactor_int[i] < magFactor_int[i-1])
                {
                    cout << "magsInput must be strictly non-decreasing" << endl;
                    MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Abort(MPI_COMM_WORLD, -1);
                }
            }
        }

        std::sort(magFactor_int.begin(), magFactor_int.end());                  // sort
        auto last = std::unique(magFactor_int.begin(), magFactor_int.end());    // get unique elements
        magFactor_int.erase(last, magFactor_int.end());

        nUniqueGroups = magFactor_int.size();

        groupRef.resize(nUniqueGroups);
        varNames.resize(nUniqueGroups);
        numVars = nUniqueGroups;

        // get indices within dataset for variables in each group
        for (int i = 0; i < setData->numVars; ++i)
        {
            for (int j = 0; j < nUniqueGroups; ++j)
            {
                if (stoi(magsToken[i]) == magFactor_int[j])
                {
                    groupRef[j].push_back(i);
                    varNames[j] += setData->varName[i];
                }
            }
        }

        for (int j = 0; j < nUniqueGroups; ++j)
        {
            if (groupRef[j].size() > 1)
            {
                varNames[j] += "Mag";
            }
        }

        if (numVars == setData->numVars)
        {
            cout << "Requested magnitude calculations, but no matching variables" << endl;
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
    }
    else
    {
        numVars = setData->numVars;
        varNames = setData->varName;
    }

    int numDOF = numVars * setData->nCells;
    pMat* dataTruthCalc;
    pMat* dataCompCalc;
    if (mags)
    {
        dataTruthCalc = new pMat(numDOF, setData->nSets, dataTruth->pG, false);
        calc_magnitude_data(dataTruth, dataTruthCalc, setData, groupRef);
        dataCompCalc = new pMat(numDOF, setData->nSets, dataComp->pG, false);
        calc_magnitude_data(dataComp, dataCompCalc, setData, groupRef);
    }
    else
    {
        dataTruthCalc = dataTruth;
        dataCompCalc = dataComp;
    }

    pMat* onesRow = new pMat(1, setData->nCells, dataTruthCalc->pG, 0, 0, 1.0, false);
    pMat* errField = new pMat(numDOF, setData->nSets, dataTruthCalc->pG, false);
    pMat* errVar = new pMat(numVars, setData->nSets, dataTruthCalc->pG, false);
    pMat* errVarNorm = new pMat(numVars, setData->nSets, dataTruthCalc->pG, false);
    pMat* norm = new pMat(numVars, setData->nSets, dataTruthCalc->pG, false);

    pMat* onesRowVars = new pMat(1, numVars, dataTruthCalc->pG, 0, 0, 1.0, false);
    pMat* errAvgNorm = new pMat(1, setData->nSets, dataTruthCalc->pG, false);
    pMat* errAvgNormP0 = new pMat(1, setData->nSets, dataTruthCalc->pG, 0, 2, 0.0, false);

    // # absolute error #
    // full field error snapshots [abs(qTruth - qComp)]
    cout << "Calculating absolute error" << endl;
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(dataTruthCalc->dataD[i] - dataCompCalc->dataD[i]);
    if (outErrField)
        setData->batchWrite(errField, errDir, "abs_err_");

    // variable error history [sum(abs(qTruth - qComp)) / nCells]
    for (int i = 0; i < numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= setData->nCells;
    errVar->write_bin(errDir + "/abs_avg_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, varNames, errDir + "/abs_avg_sum_err" + errSuffix + ".dat");

    // relative variable error history [sum(abs(qTruth - qComp)) / sum(abs(qTruth))]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(dataTruthCalc->dataD[i]);
    for (int i = 0; i < numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVarNorm->dataD[i] = errVar->dataD[i] / (norm->dataD[i] / setData->nCells);  // need to divide by nCells to negate nCells in numerator
    errVarNorm->write_bin(errDir + "/abs_rel_err" + errSuffix + ".bin");
    calc_integrated_error(errVarNorm, varNames, errDir + "/abs_rel_sum_err" + errSuffix + ".dat");

    // # L2 error #
    cout << "Calculating L2 error" << endl;
    for (int i = 0; i < errField->dataD.size(); ++i) {
        errField->dataD[i] = dataTruthCalc->dataD[i] - dataCompCalc->dataD[i];
        errField->dataD[i] = errField->dataD[i] * errField->dataD[i];
    }

    // variable error history [||qTruth - qComp|| / nCells]
    for (int i = 0; i < numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] = sqrt(errVar->dataD[i]) / setData->nCells;
    errVar->write_bin(errDir + "/l2_avg_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, varNames, errDir + "/l2_avg_sum_err" + errSuffix + ".dat");

    // relative variable error history [||qTruth - qComp|| / ||qTruth||]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = dataTruthCalc->dataD[i] * dataTruthCalc->dataD[i];
    for (int i = 0; i < numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < norm->dataD.size(); ++i)
        norm->dataD[i] = sqrt(norm->dataD[i]);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVarNorm->dataD[i] = errVar->dataD[i] / (norm->dataD[i] / setData->nCells);  // need to divide by nCells to negate nCells in numerator
    errVarNorm->write_bin(errDir + "/l2_rel_err" + errSuffix + ".bin");
    calc_integrated_error(errVarNorm, varNames, errDir + "/l2_rel_sum_err" + errSuffix + ".dat");

    // print to STDOUT
    errAvgNorm->matrix_Product('N', 'N', 1, setData->nSets, numVars, onesRowVars, 0, 0, errVarNorm, 0, 0, 1.0 / numVars, 0.0, 0, 0);
    errAvgNormP0->changeContext(errAvgNorm, false);

    int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
        for (int i = 0; i < errAvgNormP0->dataD.size(); i++) {
            cout << "Sample " << (i + 1) << endl;
            cout << "Average normalized error: " << setprecision(numeric_limits<double>::digits10) << scientific << errAvgNormP0->dataD[i] << endl;
        }
    }

    destroyPMat(onesRow, false);
    destroyPMat(errField, false);
    destroyPMat(errVar, false);
    destroyPMat(errVarNorm, false);
    destroyPMat(errAvgNorm, false);
    destroyPMat(errAvgNormP0, false);
    destroyPMat(norm, false);

}
