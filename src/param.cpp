#include "param.hpp"

paramMap::paramMap(std::string file)
{
        filename = file;
        keys.clear();
        intParam.clear();
        doubleParam.clear();
}
paramMap::paramMap(std::string file,int r)
{
        isMPI=true;
        rank=r;
        filename = file;
        keys.clear();
        intParam.clear();
        doubleParam.clear();
}
bool paramMap::getParamInt(std::string parastr, int &Param)
{
        //std::cout<<"Looking for "<<parastr<<std::endl;
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
                                std::string token=line.substr(0, eq);
                                token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                if (parastr == token)
                                {
                                        //	std::cout<<"found"<<std::endl;
                                        token=line.substr(eq + 1);
                                        token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                        Param = std::stoi(token, nullptr);
                                        return true;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        if(isMPI&& (!rank))
                throw(-1);
        if(isMPI==false)
                throw(-1);
        return false;
}
bool paramMap::getParamDouble(std::string parastr, double &Param)
{
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
                                std::string token=line.substr(0, eq);
                                token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                if (parastr == token)
                                {
                                        //std::cout<<"found"<<std::endl;
                                        token=line.substr(eq + 1);
                                        token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                        Param = std::stod(token, nullptr);
                                        return true;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        if(isMPI&& (!rank))
                throw(-1);
        if(isMPI==false)
                throw(-1);
        return false;
}

bool paramMap::getParamString(std::string parastr,std::string &Param)
{
        //std::cout<<"Looking for "<<parastr<<std::endl;

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
                                std::string token=line.substr(0, eq);
                                token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                if (parastr == token)
                                {
                                        //std::cout<<"found"<<std::endl;
                                        token=line.substr(eq + 1);
                                        token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                        Param = token;
                                        return true;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        if(isMPI&& (!rank))
                throw(-1);
        if(isMPI==false)
                throw(-1);

        return false;
}



//std::string paramMap::buildString