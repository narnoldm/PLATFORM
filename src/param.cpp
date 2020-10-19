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
                                std::string token=line.substr(0, eq);
                                token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                if (parastr == token)
                                {
                                        //	std::cout<<"found"<<std::endl;
                                        token=line.substr(eq + 1);
                                        token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                        Param = std::stoi(token, nullptr);
                                        return Param;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        if(isMPI&& (!rank))
                throw(-1);
        if(isMPI==false)
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
                                std::string token=line.substr(0, eq);
                                token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                if (parastr == token)
                                {
                                        //std::cout<<"found"<<std::endl;
                                        token=line.substr(eq + 1);
                                        token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                        Param = std::stod(token, nullptr);
                                        return Param;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        if(isMPI&& (!rank))
                throw(-1);
        if(isMPI==false)
                throw(-1);
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
                                std::string token=line.substr(0, eq);
                                token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                if (parastr == token)
                                {
                                        //std::cout<<"found"<<std::endl;
                                        token=line.substr(eq + 1);
                                        token.erase(remove_if(token.begin(),token.end(),isspace),token.end());
                                        Param = token;
                                        return Param;
                                }
                        }
                }
        }
        std::cout << "param not found " << parastr << std::endl;
        if(isMPI&& (!rank))
                throw(-1);
        if(isMPI==false)
                throw(-1);
}



//std::string paramMap::buildString