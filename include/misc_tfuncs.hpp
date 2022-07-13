
#ifndef MISC_TFUNC_H
#define MISC_TFUNC_H

template <class T>
void printASCIIVecP0(std::string fname, std::vector<T> &Mat, int N)
{
        FILE *fid;
        fid = fopen(fname.c_str(), "w");
        fprintf(fid, "%d\n", N);
        for (int i = 0; i < N; i++)
        {
                if (typeid(T) == typeid(double))
                        fprintf(fid, "%.9E\n", Mat[i]);
                if (typeid(T) == typeid(int))
                        fprintf(fid, "%d\n", (int)Mat[i]);
        }
        fclose(fid);
}

#endif