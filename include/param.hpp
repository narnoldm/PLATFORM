
#ifndef PARAM_H
#define PARAM_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include <algorithm>

///
/// paramMap is the class used for reading in the input file parameters
/// This class doesn't strictly need MPI as each processor will read the input files independently
/// which can be helpful when used within logic blocks
///
class paramMap
{
public:
	std::string filename;	 ///< The paramMap filename to read from
	int rank;				 ///< The local rank of the processor
	bool isMPI = false;		 ///< Whether or not this is the 0th rank
	bool catchMisInp = true; ///< Whether or not to catch missing input (no implemented)

public:
	paramMap(std::string file);																///< Constructor
	paramMap(std::string file, int r);														///< Constructor with MPI rank
	bool getParamInt(std::string key, int &Param);											///< Get an integer parameter
	bool getParamInt(std::string key, int &Param, const int defaultVal);					///< Get an integer parameter with a default value
	bool getParamDouble(std::string key, double &Param);									///< Get a double parameter
	bool getParamDouble(std::string key, double &Param, const double defaultVal);			///< Get a double parameter with a default value
	bool getParamString(std::string key, std::string &Param);								///< Get a string parameter
	bool getParamString(std::string key, std::string &Param, const std::string defaultVal); ///< Get a string parameter with a default value
	bool getParamBool(std::string key, bool &Param);										///< Get a bool parameter
	bool getParamBool(std::string key, bool &Param, const bool defaultVal);					///< Get a bool parameter with a default value

	template <class T>
	bool parseSep(const std::string &input, std::string sep, std::vector<T> &tokens); ///< Parse a string into a vector of tokens using seperator
};

#endif