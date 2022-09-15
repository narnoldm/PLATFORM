
\mainpage PLATFORM
 



\section intro_sec Introduction
 
In the scientific model reduction community a variety of 
decomposition techniques are used for analysis of complex
spatio-temporal systems as well as reduced-order model 
development. Application of these methods to more complex
problems require memory usage that often exceeds that of 
single machines. Parallel Linear Algebra Tool FOr Reduced Modeling  (PLATFORM) is a MPI 
C++ code built upon the ScaLAPACK toolset of distributed
Linear Algebra to address these challenges.



\section install_sec Installation Instructions

PLATFORM uses CMake to configure and build. Detailed instructions can be found [here](@ref installation).


\section example Example Cases
 
 An example case is provide in the examples directory for how to use the code to read in binary data and perform a POD decomposition [Example 1](@ref example1)



 
