#include "error_funcs.hpp"

using namespace :: std;

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

void calc_abs_and_l2_error(pMat* dataTruth, pMat* dataComp, tecIO* setData, string errDir, string errSuffix, bool outErrField) {

    pMat* onesRow = new pMat(1, setData->nCells, dataTruth->pG, 0, 0, 1.0, false);
    pMat* errField = new pMat(setData->nPoints, setData->nSets, dataTruth->pG, false);
    pMat* errVar = new pMat(setData->numVars, setData->nSets, dataTruth->pG, false);
    pMat* errVarNorm = new pMat(setData->numVars, setData->nSets, dataTruth->pG, false);
    pMat* norm = new pMat(setData->numVars, setData->nSets, dataTruth->pG, false);

    pMat* onesRowVars = new pMat(1, setData->numVars, dataTruth->pG, 0, 0, 1.0, false);
    pMat* errAvgNorm = new pMat(1, setData->nSets, dataTruth->pG, false);
    pMat* errAvgNormP0 = new pMat(1, setData->nSets, dataTruth->pG, 0, 2, 0.0, false);

    // # absolute error #
    // full field error snapshots [abs(qTruth - qComp)]
    cout << "Calculating absolute error" << endl;
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(dataTruth->dataD[i] - dataComp->dataD[i]);
    if (outErrField)
        setData->batchWrite(errField, errDir, "abs_err_");

    // variable error history [sum(abs(qTruth - qComp)) / nCells]
    for (int i = 0; i < setData->numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] /= setData->nCells;
    errVar->write_bin(errDir + "/abs_avg_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, setData->varName, errDir + "/abs_avg_sum_err" + errSuffix + ".dat");

    // relative variable error history [sum(abs(qTruth - qComp)) / sum(abs(qTruth))]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = abs(dataTruth->dataD[i]);
    for (int i = 0; i < setData->numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVarNorm->dataD[i] = errVar->dataD[i] / (norm->dataD[i] / setData->nCells);  // need to divide by nCells to negate nCells in numerator
    errVarNorm->write_bin(errDir + "/abs_rel_err" + errSuffix + ".bin");
    calc_integrated_error(errVarNorm, setData->varName, errDir + "/abs_rel_sum_err" + errSuffix + ".dat");

    // # L2 error #
    cout << "Calculating L2 error" << endl;
    for (int i = 0; i < errField->dataD.size(); ++i) {
        errField->dataD[i] = dataTruth->dataD[i] - dataComp->dataD[i];
        errField->dataD[i] = errField->dataD[i] * errField->dataD[i];
    }

    // variable error history [||qTruth - qComp|| / nCells]
    for (int i = 0; i < setData->numVars; ++i)
        errVar->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVar->dataD[i] = sqrt(errVar->dataD[i]) / setData->nCells;
    errVar->write_bin(errDir + "/l2_avg_err" + errSuffix + ".bin");
    calc_integrated_error(errVar, setData->varName, errDir + "/l2_avg_sum_err" + errSuffix + ".dat");

    // relative variable error history [||qTruth - qComp|| / ||qTruth||]
    // re-use errField to save memory
    for (int i = 0; i < errField->dataD.size(); ++i)
        errField->dataD[i] = dataTruth->dataD[i] * dataTruth->dataD[i];
    for (int i = 0; i < setData->numVars; ++i)
        norm->matrix_Product('N', 'N', 1, setData->nSets, setData->nCells, onesRow, 0, 0, errField, i * setData->nCells, 0, 1.0, 0.0, i, 0);
    for (int i = 0; i < norm->dataD.size(); ++i)
        norm->dataD[i] = sqrt(norm->dataD[i]);
    for (int i = 0; i < errVar->dataD.size(); ++i)
        errVarNorm->dataD[i] = errVar->dataD[i] / (norm->dataD[i] / setData->nCells);  // need to divide by nCells to negate nCells in numerator
    errVarNorm->write_bin(errDir + "/l2_rel_err" + errSuffix + ".bin");
    calc_integrated_error(errVarNorm, setData->varName, errDir + "/l2_rel_sum_err" + errSuffix + ".dat");

    // print to STDOUT
    errAvgNorm->matrix_Product('N', 'N', 1, setData->nSets, setData->numVars, onesRowVars, 0, 0, errVarNorm, 0, 0, 1.0 / setData->numVars, 0.0, 0, 0);
    errAvgNormP0->changeContext(errAvgNorm);

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
