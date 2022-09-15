Installation                 {#installation}
=============



PLATFORM is a C++ code and is configured and built using the CMake system. 


Prerequisites
-------------

The Following are required to build PLATFORM:

1. BLAS libraries (e.g. OpenBLAS, MKL) 
2. LAPACK libraries 
3. MPI libraries (e.g. OpenMPI, MPICH)


The PLATFORM CMakeLists will attempt to download and build the following libraries(web access is required):

1. Boost (used for the TecIO libraries)
2. Eigen (used for some single process operations)
3. ScaLAPACK (Only if not provided by Math Kernal Library) 


Building PLATFORM
-----------------

To launch the cmake configuration and build process, create a build directory and run cmake from within the build directory.
~~~~
mkdir build
cd build
cmake .. 
~~~~

This should complete without issue. If not, please see the Common Issues section below.

PLATFORM is usually built on a cluster where the compiler is set in the environment via modules. 
If you are compiling on a local machine with multiple compilers you may need to explicitly link the compiler using

Example for intel MPI compiler
~~~
cmake .. -DCMAKE_CXX_COMPILER=mpiicpc -DCMAKE_C_COMPILER=mpiicc -DCMAKE_Fortran_COMPILER=mpiifort
~~~

After the configuration concludes you can build the code using 

~~~
make 
~~~
or in parallel
~~~
make -j 4
~~~

It is advised to use a serial build to make parsing any errors easier. 


Common Issues
-------------