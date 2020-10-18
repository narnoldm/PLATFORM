#include "param.hpp"

paramMap::paramMap(std::string file)
{
        filename = file;
        keys.clear();
        intParam.clear();
        doubleParam.clear();
}

int paramMap::getParamInt(std::string parastr)
{
        //std::cout<<"Looking for "<<parastr<<std::endl;
        int Param;
        std::fstream inFile;
        inFile.open(filename);
        std::string line;

        //std::stringstream token;

        while (!inFile.eof())
        {
                std::getline(inFile, line);
                //std::cout<<line<<std::endl;
                if (line.length() != 0)
                {
                        //find =
                        int eq = line.find("=");
                        //std::cout<<eq<<std::endl;
                        if (line.length() != eq)
                        {
                                //std::cout<<line.substr(0,eq)<<std::endl;
                                //token>>line.substr(0,eq);
                                if (parastr == line.substr(0, eq))
                                {
                                        //	std::cout<<"found"<<std::endl;
                                        Param = std::stoi(line.substr(eq + 1), nullptr);
                                        return Param;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        throw(-1);
        return 0;
}
double paramMap::getParamDouble(std::string parastr)
{
        //std::cout<<"Looking for "<<parastr<<std::endl;
        double Param;
        std::fstream inFile;
        inFile.open(filename);
        std::string line;

        //std::stringstream token;

        while (!inFile.eof())
        {
                std::getline(inFile, line);
                //std::cout<<line<<std::endl;
                if (line.length() != 0)
                {
                        //find =
                        int eq = line.find("=");
                        //std::cout<<eq<<std::endl;
                        if (line.length() != eq)
                        {
                                //std::cout<<line.substr(0,eq)<<std::endl;
                                //token>>line.substr(0,eq);
                                if (parastr == line.substr(0, eq))
                                {
                                        //std::cout<<"found"<<std::endl;
                                        Param = std::stod(line.substr(eq + 1), nullptr);
                                        return Param;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        throw(-1);
        return 0;
}

std::string paramMap::getParamString(std::string parastr)
{
        //std::cout<<"Looking for "<<parastr<<std::endl;
        std::string Param;
        std::fstream inFile;
        inFile.open(filename);
        std::string line;

        //std::stringstream token;

        while (!inFile.eof())
        {
                std::getline(inFile, line);
                //std::cout<<line<<std::endl;
                if (line.length() != 0)
                {
                        //find =
                        int eq = line.find("=");
                        //std::cout<<eq<<std::endl;
                        if (line.length() != eq)
                        {
                                //std::cout<<line.substr(0,eq)<<std::endl;
                                //token>>line.substr(0,eq);
                                if (parastr == line.substr(0, eq))
                                {
                                        //std::cout<<"found"<<std::endl;
                                        Param = line.substr(eq + 1);
                                        return Param;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        throw(-1);
        return " ";
}



//std::string paramMap::buildString