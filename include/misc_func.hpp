
#ifndef MISC_FUNC_H
#define MISC_FUNC_H

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <type_traits>
#include <typeinfo>
#include <assert.h>
#include <algorithm>
#include <cctype>




void printASCIIVecP0(std::string fname, double *Mat, int N);

bool to_bool(std::string str);

void tokenparse(const std::string &input, std::string sep, std::vector<std::string> &tokens);




#endif