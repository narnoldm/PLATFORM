cmake_minimum_required(VERSION 3.16.2)
project(PDP VERSION 3.0)

include(ExternalProject)
find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 17)
set(TECDIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_Party/tecio/teciosrc)


set(EIGEN_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external/eigen)
ExternalProject_Add(Eigen3
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} 
    UPDATE_COMMAND ""
    CONFIGURE_COMMAND ""
    CMAKE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
)
include_directories(${CMAKE_BINARY_DIR}/Eigen3-prefix/src/Eigen3)


set(GTEST_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external/googleTest)
ExternalProject_Add(googletest
    GIT_REPOSITORY https://github.com/google/googletest
    CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_INSTALL_PREFIX=${GTEST_INSTALL_LOCATION}
)
include_directories(${GTEST_INSTALL_LOCATION}/include)
link_directories(${GTEST_INSTALL_LOCATION}/lib)








if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
message("${CMAKE_CXX_COMPILER_ID} detected")
set(isINTEL 1)
message("MKLROOT $ENV{MKLROOT}")
set(MKLROOT $ENV{MKLROOT})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -DUSE_MKL")
endif()




add_subdirectory(src)
add_subdirectory(CASLAB_execs)








enable_testing()
add_subdirectory(tst)



