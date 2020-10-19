
#ifndef PARAM_H
#define PARAM_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <mpi.h>
#include <algorithm>


class paramMap
{
  public:
	std::string filename;
	std::vector<std::string> keys;
	std::vector<int> intParam;
	std::vector<double> doubleParam;
	std::vector<std::string> stringParam;
	int rank;
	bool isMPI=false;

  public:
	paramMap(std::string file);
	paramMap(std::string file,int r);
	int getParamInt(std::string key);
	double getParamDouble(std::string key);
	std::string getParamString(std::string key);
};


#endif