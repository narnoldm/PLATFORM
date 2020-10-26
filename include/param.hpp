
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
	bool catchMisInp=true;

  public:
	paramMap(std::string file);
	paramMap(std::string file,int r);
	bool getParamInt(std::string key,int & Param);
	bool getParamDouble(std::string key,double & Param);
	bool getParamString(std::string key,std::string & Param);
	bool getParamBool(std::string key, bool &Param);

	template <class T>
	bool parseSep(const std::string &input, std::string sep, std::vector<T> &tokens);

};


#endif