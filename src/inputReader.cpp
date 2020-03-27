#include "inputReader.hpp"



inputReader::inputReader(string file)
{
        cout << "file is :" << file << endl;
        ifile = file;
        ScanFile();
}
inputReader::~inputReader()
{
        ifile.clear();
}
bool inputReader::ScanFile()
{
        fstream inFile;
        inFile.open(ifile.c_str(), std::fstream::in);
        string line;
        assert(inFile.good());

        enum blocks
        {
                info,
                input,
                operation,
                output
        };
        int cBlock = info;
        int comment;
        while (!inFile.eof())
        {
                getline(inFile, line);
                comment = line.find("!");
                if (comment != string::npos)
                        line = line.substr(0, comment);
                if (line.empty())
                {
                        continue;
                }
                if (line == "#INFO")
                {
                        cout << line << endl;
                        continue;
                }
                if (line == "#INPUT")
                {
                        cout << line << endl;
                        cBlock = input;
                        continue;
                }
                if (line == "#OPERATIONS")
                {
                        cout << line << endl;
                        cBlock = operation;
                        continue;
                }
                if (line == "#OUTPUT")
                {
                        cout << line << endl;
                        cBlock = output;
                        continue;
                }
                switch (cBlock)
                {
                case info:
                        sys.ScanInfo(line);
                        break;
                case input:
                        inp.ScanInput(line);
                        break;
                case operation:
                        op.ScanOperations(line);
                        break;
                case output:
                        out.ScanInput(line);
                        break;
                }
        }
        inFile.close();
        cout << sys;
        cout << inp;
        cout << op;
        op.assignMats(inp);
        cout << op;
        cout << out;

        int i = 0;
        if(!inp.checkMats())
        {
                cout << "not all mats defined" << endl;
                inp.assignInputInfo();
                op.inferDim();

                i++;
                if (i > 10)
                {
                        cout << "matricies undefined after 10 iterations: check input file" << endl;
                        throw(-1);
                }
        }

        if(out.matList.size()==0)
        {
                cout<<"No outputs defined. Add an output"<<endl;
                throw(-1);
        }

        cout << "Yay all orperations and matracies defined" << endl;
}