
#ifndef MISC_TFUNC_H
#define MISC_TFUNC_H

template <class T>
void printASCIIVecP0(std::string fname, std::vector<T> &Mat, int N)
{
        FILE *fid;
        fid = fopen(fname.c_str(), "w");
        int ierr = fprintf(fid, "%d\n", N);
        for (int i = 0; i < N; i++)
        {
                if (typeid(T) == typeid(double))
                        ierr = fprintf(fid, "%.9E\n", (double)Mat[i]);
                if (typeid(T) == typeid(int))
                        ierr = fprintf(fid, "%d\n", (int)Mat[i]);
        }
        fclose(fid);
}

#endif