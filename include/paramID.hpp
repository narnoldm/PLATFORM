

#ifndef PARAMID_H
#define PARAMID_H


#include <iostream>
#include <string>
#include <typeinfo>
#include "misc_func.hpp"

using namespace::std;


template <class T>
class paramID
{
public:
    string key;
    bool set;
    T value;

    paramID();
    paramID(string kIn, T vIn);
    paramID(string kIn);
    ~paramID();

    bool checkInfo(string line);
};

template <typename T>
std::ostream &operator<<(std::ostream &os, const paramID<T> &pID);


template <typename T>
ostream &operator<<(std::ostream &os, const paramID<T> &pID)
{
        cout << pID.key << " = " << pID.value << " and is of type: " << typeid(pID.value).name() << endl;
        return os;
}
template <class T>
paramID<T>::paramID()
{
        key = "";
        value = 0;
        set = false;
};
template <class T>
paramID<T>::paramID(string kIn, T vIn)
{
        key = kIn;
        value = vIn;
        set = true;
};

template <class T>
paramID<T>::paramID(string kIn)
{
        key = kIn;
        value = 0;
        set = false;
};

template <class T>
paramID<T>::~paramID(){};

template <class T>
bool paramID<T>::checkInfo(string line)
{
        int eq = line.find("=");
        if (line.substr(0, eq) != key)
                return 0;
        else
        {
                if (typeid(T) == typeid(int))
                        value = stoi(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(double))
                        value = stod(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(float))
                        value = stof(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(long))
                        value = stol(line.substr(eq + 1), nullptr);
                else if (typeid(T) == typeid(bool))
                        value = to_bool(line.substr(eq + 1));
                else
                {
                        cout << "string to type not defined for " << typeid(T).name() << endl;
                        return 0;
                }
                set = 1;
                return 1;
        }
};



#endif