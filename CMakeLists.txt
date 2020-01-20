cmake_minimum_required(VERSION 3.16.2)
project(PDP VERSION 3.0)

message("MKLROOT $ENV{MKLROOT}")
set(MKLROOT $ENV{MKLROOT})
set(CMAKE_CXX_STANDARD 11)


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
message("${CMAKE_CXX_COMPILER_ID} detected")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()



add_library(pMat STATIC src/pMat.cpp)
target_include_directories(pMat PRIVATE include)
add_library(parser STATIC src/parser.cpp)
target_include_directories(parser PRIVATE include)

target_link_libraries(parser pMat)

target_link_directories(pMat PUBLIC ${MKLROOT}/lib/intel64/)
#target_link_libraries(pMat mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64) 
target_link_libraries(pMat ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 

target_include_directories(pMat PRIVATE include ${MKLROOT}/include)


add_executable(parser_drive src/parser_drive.cpp)
target_include_directories(parser_drive PRIVATE include)
target_link_libraries(parser_drive parser)



enable_testing()

add_test(NAME test1_SVD_bin_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/parser_drive)
add_test(NAME test1_SVD_bin_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/parser_drive)












