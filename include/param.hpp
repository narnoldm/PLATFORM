#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

class paramMap
{
  public:
	std::string filename;
	std::vector<std::string> keys;
	std::vector<int> intParam;
	std::vector<double> doubleParam;
	std::vector<std::string> stringParam;

  public:
	paramMap(std::string file);
	int getParamInt(std::string key);
	double getParamDouble(std::string key);
	std::string getParamString(std::string key);
};
