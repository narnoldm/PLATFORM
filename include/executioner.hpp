#ifndef EXECUTIONER_H
#define EXECUTIONER_H



#include "inputReader.hpp"
#include<vector>
#include<iostream>


class executioner
{
    public: 
    inputReader * inpFile;
    vector<PGrid*> pGs;


    executioner();
    executioner(inputReader * iF);
    ~executioner();

    void init();
    void exec_all();
    void exec(int opID);
    void output();
    void clear();
    void create_matricies();
};


#endif