
#ifndef MISC_TFUNC_H
#define MISC_TFUNC_H

template <class T>
void printASCIIVecP0(std::string fname, std::vector<T> & Mat, int N)
{
        FILE *fid;
        fid = fopen(fname.c_str(), "w");
        for (int i = 0; i < N; i++)
        {
                if(typeid(T)==typeid(double))
                        fprintf(fid, "%d %.9E\n", i, Mat[i]);
                if(typeid(T)==typeid(int))
                        fprintf(fid, "%d %d\n", i, Mat[i]);
        }
        fclose(fid);
}


#endif