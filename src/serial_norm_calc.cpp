#include "inputReader.hpp"

#include <limits>
typedef numeric_limits<double> dbl;

void prepConsVarProc(tecIO*, int&, int&, int&, int&, vector<int>&, vector<int>&, vector<vector<int>>&, vector<int>&);
void calcConsVars(double, double, double, vector<double>, double, double, vector<double>&);

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); //Basic MPI intialization
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    //silence output for non root MPI processes
    std::ofstream sink("/dev/null");
    streambuf *strm_buffer = cout.rdbuf();
    int debug_proc = 0;
    int calcConsv = 0;
    string avgFile;
    if(argc > 2)
    {
        debug_proc = atoi(argv[2]);
        if (argc > 3) {
            avgFile = argv[3];
            
            // calculate conservative normalization values
            if (argc > 4) {
                calcConsv = atoi(argv[4]);
            }
        }

    }
    if (rank != debug_proc)
    {
        std::cout.rdbuf(sink.rdbuf());
    }
    string input = argv[1];
    cout << "input string is: " << input << endl;
    vector<string> tokens;
    tokenparse(input, "|", tokens);

    tecIO *dataset1 = new tecIO(tokens);
    //check for non standard input file name

    vector<double> data,average,norm,sum;
    int N = dataset1->nCells*dataset1->numVars;
    int M = dataset1->nSets;
    data.resize(N);
    average.resize(N);

    cout.precision(dbl::max_digits10);

    if (argc < 4) {
        double val = 0.0;
        sum.resize(N);

        for(int i = dataset1->snap0; i <= dataset1->snapF; i = i + dataset1->snapSkip)
        {
            dataset1->readSingle(i, data.data());
            
            for(int j = 0; j < data.size(); j++)
                average[j] += data[j];
        }
        for(int j = 0; j < data.size(); j++)
            average[j] /= M;

        cout << average[0] << endl;

        for(int i = dataset1->snap0; i <= dataset1->snapF; i = i + dataset1->snapSkip)
        {
            dataset1->readSingle(i,data.data());
            for(int j = 0; j < data.size(); j++)
                sum[j] += (data[j] - average[j])*(data[j] - average[j]);
        }
        cout << sum[0] <<endl;
        for(int j = 0; j < data.size(); j++)
                sum[j] /= M;

        cout << sum[0] <<endl;
        for(int j = 0; j < data.size(); j++)
            val += sum[j];

	    cout << val << endl;
        val /= N;
        val = sqrt(val);

        cout << val << endl;

    } else {
        // load average file
        dataset1->readAvg(avgFile);

        vector<int> normFactor_int;
        int nUniqueGroups;
        vector<vector<int>> groupRef;
        int pIdx, tIdx, densityIdx, enthalpyIdx;
        vector<int> vIdxs;
        vector<int> scalarIdxs;

        // TODO: there should be a way to account for preset (positive) norm factors.
        // could theoretically just ignore those, as this implies we already know the norm factor to use.
        //      probably just exclude from groupRef. May not be worth time, just let calcs roll


        // prep group references for primitive variables manually
        // set up references for calculating conserved variables
        if (calcConsv) {

            prepConsVarProc(dataset1, pIdx, tIdx, densityIdx, enthalpyIdx, vIdxs, scalarIdxs, groupRef, normFactor_int);
            nUniqueGroups = normFactor_int.size();

        // otherwise do this in a more general manner
        } else {
            // get sorted list of unique norm factors
            normFactor_int.resize((dataset1->numVars));
            std::copy(dataset1->normFactor.begin(), dataset1->normFactor.end(), normFactor_int.begin());    // convert to integers
            std::sort(normFactor_int.begin(),normFactor_int.end());                                   // sort
            auto last = std::unique(normFactor_int.begin(), normFactor_int.end());                    // get unique elements
            normFactor_int.erase(last, normFactor_int.end());

            nUniqueGroups = normFactor_int.size();

            groupRef.resize(nUniqueGroups);

            // get indices within dataset for variables in each group
            for (int i = 0; i < dataset1->numVars; ++i) {
                for (int j = 0; j < nUniqueGroups; ++j) {
                    if (int(dataset1->normFactor[i]) == normFactor_int[j]) {
                        groupRef[j].push_back(i);
                    }
                }
            }
        }

        cout << "number unique groups: " << nUniqueGroups << endl;

        int dataSizeGroup = dataset1->nCells*nUniqueGroups;
        vector<double> dataAvg(dataSizeGroup);
        vector<double> dataSnap(dataSizeGroup);
        vector<double> sum_consv;
        sum.resize(dataSizeGroup);

        // compute magnitude data for average file
        int idx_full, idx_group;
        double magVal;
        cout << "Calculating magnitudes for average file..." << endl;
        for (int j = 0; j < nUniqueGroups; ++j) {

            // compute magnitude for each group
            for (int iCell = 0; iCell < dataset1->nCells; ++iCell) {

                magVal = 0;
                // don't calculate magnitude if it's a single variable
                if (groupRef[j].size() == 1) {
                    magVal = dataset1->average[groupRef[j][0]*dataset1->nCells + iCell];
                } else {
                    for (int k = 0; k < groupRef[j].size(); ++k) {
                        idx_full = groupRef[j][k]*dataset1->nCells + iCell;
                        magVal += (dataset1->average[idx_full])*(dataset1->average[idx_full]);
                    }
                    magVal = sqrt(magVal);
                }
                dataAvg[j*dataset1->nCells + iCell] = magVal;
            }
        }

        // compute conserved variables for average
        // probably takes a speed-hit from non-contiguous memory reads, but no two ways about it,
        //          need to get all primitive vars and dens/enth from a given cell
        // TODO: shunt this off to a function to make it less confusing/cluttered
        double pressure, velMag, temperature, density, enthalpy;
        vector<double> scalars;
        vector<double> consvVals;
        vector<double> dataAvg_consv;
        if (calcConsv) {
            scalars.resize(scalarIdxs.size());
            consvVals.resize(nUniqueGroups);
            dataAvg_consv.resize(dataSizeGroup);
            for (int iCell = 0; iCell < dataset1->nCells; ++iCell) {
                // collect data necessary to compute conservative variables
                pressure  = dataAvg[iCell];
                velMag = dataAvg[dataset1->nCells + iCell];
                temperature   = dataAvg[2*dataset1->nCells + iCell];
                for (int i = 0; i < scalarIdxs.size(); ++i) {
                    scalars[i] = dataAvg[(3+i)*dataset1->nCells + iCell];
                }
                density  = dataset1->average[densityIdx*dataset1->nCells + iCell];
                enthalpy = dataset1->average[enthalpyIdx*dataset1->nCells + iCell];

                // calculate conserved variables, distribute into vector
                calcConsVars(pressure, temperature, velMag, scalars, density, enthalpy, consvVals);
                for (int i = 0; i < nUniqueGroups; ++i) {
                    dataAvg_consv[i*dataset1->nCells + iCell] = consvVals[i];
                }
            }
        }

        vector<double> dataSnap_consv;
        if (calcConsv) {
            dataSnap_consv.resize(dataSizeGroup);
            sum_consv.resize(dataSizeGroup);
        }

        // read through each snapshot
        cout << "Calculating L2 norm for mean-subtracted values..." << endl;
        for(int i = dataset1->snap0; i <= dataset1->snapF; i= i + dataset1->snapSkip) {
            dataset1->readSingle(i,data.data());

            // compute running sum for each group
            for (int j = 0; j < nUniqueGroups; ++j) {

                // compute magnitude for each group
                for (int iCell = 0; iCell < dataset1->nCells; ++iCell) {

                    magVal = 0;
                    if (groupRef[j].size() == 1) {
                        magVal = data[groupRef[j][0]*dataset1->nCells + iCell];
                    } else {
                        for (int k = 0; k < groupRef[j].size(); ++k) {
                            idx_full = groupRef[j][k]*dataset1->nCells + iCell;
                            magVal += (data[idx_full])*(data[idx_full]);
                        }
                        magVal = sqrt(magVal);
                    }
                    idx_group = j*dataset1->nCells + iCell;
                    dataSnap[idx_group] = magVal;
                    sum[idx_group] += (magVal - dataAvg[idx_group])*(magVal - dataAvg[idx_group]);

                }
            }

            // compute running sum for conservative variable magnitudes
            if (calcConsv) {
                for (int iCell = 0; iCell < dataset1->nCells; ++iCell) {
                    // collect data necessary to compute conservative variables
                    pressure  = dataSnap[iCell];
                    velMag = dataSnap[dataset1->nCells + iCell];
                    temperature   = dataSnap[2*dataset1->nCells + iCell];
                    for (int i = 0; i < scalarIdxs.size(); ++i) {
                        scalars[i] = dataSnap[(3+i)*dataset1->nCells + iCell];
                    }
                    density  = data[densityIdx*dataset1->nCells + iCell];
                    enthalpy = data[enthalpyIdx*dataset1->nCells + iCell];

                    // calculate conserved variables, compute square, distribute into vector
                    calcConsVars(pressure, temperature, velMag, scalars, density, enthalpy, consvVals);
                    for (int i = 0; i < nUniqueGroups; ++i) {
                        idx_group = i*dataset1->nCells + iCell;
                        sum_consv[idx_group] += (consvVals[i] - dataAvg_consv[idx_group])*(consvVals[i] - dataAvg_consv[idx_group]);
                    }
                }
            }
         }

        // compute final normalization factors
        vector<double> val(nUniqueGroups);
        for (int i = 0; i < nUniqueGroups; ++i) {
            for (int j = 0; j < dataset1->nCells; ++j) {
                val[i] += sum[i*dataset1->nCells + j];
            }
            val[i] /= M*(dataset1->nCells);
            val[i] = sqrt(val[i]);
        }

        // print in a user-friendly fashion
        for (int i = 0; i < dataset1->numVars; ++i) {

            for (int j = 0; j < nUniqueGroups; ++j) {
                if (int(dataset1->normFactor[i]) == normFactor_int[j]) {
                    cout << dataset1->varName[i] << ": " << val[j] << endl;
                    break;
                }
            }
        }
        
        // repeat for conserved variables if necessary
        vector<double> val_consv(nUniqueGroups);
        if (calcConsv) {
            for (int i = 0; i < nUniqueGroups; ++i) {
                for (int j = 0; j < dataset1->nCells; ++j) {
                    val_consv[i] += sum_consv[i*dataset1->nCells + j];
                }
                val_consv[i] /= M*(dataset1->nCells);
                val_consv[i] = sqrt(val_consv[i]);
            }

            // ordering is known
            cout << "Density: " << val_consv[0] << endl;
            cout << "Momentum: " << val_consv[1] << endl;
            cout << "Energy: " << val_consv[2] << endl;
            for (int i = 0; i < scalarIdxs.size(); ++i) {
                cout << "Density-weighted " << dataset1->varName[scalarIdxs[i]] << ": " << val_consv[3+i] << endl;
            }

        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cout.rdbuf(strm_buffer);
    MPI_Finalize();
    return 0;
}

// prep index references for computing conserved variables
// pretty much entirely hard-coded, not sure of a better way to determine pressure, temperature, etc.
//      without very explicit, tedious user input
// this is strictly to get variable ordering into a logical ordering, this function makes the cosnv var calcs  
//      totally agnostic to the order in which the user inputs the variables
// TODO: assert that all primitive variables have been included, somehow. Otherwise calcs will fail
void prepConsVarProc(tecIO* dataset, int& pressIdx, int& tempIdx, int& densIdx, int& enthIdx,
                     vector<int>& velIdxs, vector<int>& transScalarIdxs,
                     vector<vector<int>>& refIdxs, vector<int>& normFacs) {

    // get indices of variables
    for (int i = 0; i < dataset->numVars; ++i) {

        // transported scalars (ending in "_mf" or is a flamelet variable)
        if (dataset->varName[i].length() >= 3) {
            if ((dataset->varName[i].substr(dataset->varName[i].length() - 3) == "_mf") ||
                    (dataset->varName[i] == "Flamelet_Scalar_Mean") ||
                    (dataset->varName[i] == "Flamelet_Scalar_Variance") ||
                    (dataset->varName[i] == "Flamelet_Parameter")) {
                transScalarIdxs.push_back(i);
                continue;
            }
        }

        // velocity
        if ((dataset->varName[i] == "U") ||
                (dataset->varName[i] == "V") ||
                (dataset->varName[i] == "W")) {
            velIdxs.push_back(i);
            continue;
        }

        // remaining variables
        if (dataset->varName[i] == "Static_Pressure") {
            pressIdx = i;
            continue;
        }
        if (dataset->varName[i] == "Temperature") {
            tempIdx = i;
            continue;
        }
        if (dataset->varName[i] == "Density") {
            densIdx = i;
            continue;
        }
        if (dataset->varName[i] == "Enthalpy") {
            enthIdx = i;
            continue;
        }
    }

    // strictly order group IDs in the desired order,
    //      i.e. pressure, velocity, temperature, transported scalars
    vector<int> pressVec; pressVec.push_back(pressIdx);
    vector<int> tempVec; tempVec.push_back(tempIdx);
    vector<int> scalarVec(1);
    refIdxs.push_back(pressVec);
    refIdxs.push_back(velIdxs);
    refIdxs.push_back(tempVec);
    for (int i = 0; i < transScalarIdxs.size(); ++i) {
        scalarVec[0] = transScalarIdxs[i];
        refIdxs.push_back(scalarVec);
    }

    // there's definitely a more elegant way of doing this
    for (int i = 0; i < dataset->numVars; ++i) {
        if (dataset->varName[i] == "Static_Pressure") {
            normFacs.push_back(dataset->normFactor[i]);
            break;
        }
    }

    // assume that all velocity variable have same group ID
    for (int i = 0; i < dataset->numVars; ++i) {
        if (dataset->varName[i] == "U") {
            normFacs.push_back(dataset->normFactor[i]);
            break;
        }  
    }
    for (int i = 0; i < dataset->numVars; ++i) {
        if (dataset->varName[i] == "Temperature") {
            normFacs.push_back(dataset->normFactor[i]);
            break;
        }
    }

    // assumes that all transported scalars have different group IDs
    for (int i = 0; i < dataset->numVars; ++i) {
        if (dataset->varName[i].length() >= 3) {
            if ((dataset->varName[i].substr(dataset->varName[i].length() - 3) == "_mf") ||
                    (dataset->varName[i] == "Flamelet_Scalar_Mean") ||
                    (dataset->varName[i] == "Flamelet_Scalar_Variance") ||
                    (dataset->varName[i] == "Flamelet_Parameter")) {
                normFacs.push_back(dataset->normFactor[i]);
            }
        }
    }

}

// compute conservative variables from primitives, density, enthalpy
void calcConsVars(double press, double temp, double vMag, vector<double> transScalars,
                  double dens, double enth, vector<double> &consvVars) {

    // ordered by density, momentum, total enthalpy, density-weighted transported scalars
    consvVars[0] = dens;
    consvVars[1] = dens*vMag;
    consvVars[2] = dens*enth - press;
    for (int i = 0; i < transScalars.size(); ++i)
        consvVars[3+i] = dens*transScalars[i];

}