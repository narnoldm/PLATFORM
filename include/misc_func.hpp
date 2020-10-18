
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




template <class T>
void printASCIIVecP0(std::string fname, std::vector<T> & Mat, int N);
#include "misc_tfuncs.hpp"

bool to_bool(std::string str);

void tokenparse(const std::string &input, std::string sep, std::vector<std::string> &tokens);

void readMat(std::string filename, std::vector<int> & Mat);

void writeMat(std::string filename,int m,int n, std::vector<int> &Mat);

#endif