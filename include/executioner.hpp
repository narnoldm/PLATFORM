#ifndef EXECUTIONER_H
#define EXECUTIONER_H



#include "inputReader.hpp"
#include<vector>
#include<iostream>

/***
 * Executioner is the top level parser. It holds the inputReader object and the set of required process Grid buffers
 * It also ideally should be the one calling operation functions
 * */
class executioner
{
    public: 
    ///input file parser
    inputReader * inpFile;
    ///Process Grids used by all data structures
    vector<PGrid*> pGs;

    ///Default Contructor
    executioner();
    ///Contructor with passed input file
    executioner(inputReader * iF);
    ~executioner();


    ///Will attempt to allocate (and load for input type) all required data structures
    void init();
    ///Will attempt to execute each of the operations
    void exec_all();
    ///Execute the operation of given ID
    void exec(int opID);
    ///Attempt to output all output tags
    void output();
    ///Attempt to deallocate all data structures
    void clear();
    ///Attempt to create all matracies
    void create_matricies();
};


#endif