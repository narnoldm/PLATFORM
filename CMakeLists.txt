cmake_minimum_required(VERSION 3.16.2)
project(PDP VERSION 3.0)

message("MKLROOT $ENV{MKLROOT}")
set(MKLROOT $ENV{MKLROOT})
set(CMAKE_CXX_STANDARD 11)


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
message("${CMAKE_CXX_COMPILER_ID} detected")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()

include_directories(include/)


add_library(pMat STATIC src/pMat.cpp)
add_library(parser STATIC src/parser.cpp)

target_link_libraries(parser pMat)

target_link_directories(pMat PUBLIC ${MKLROOT}/lib/intel64/)
#target_link_libraries(pMat mkl_scalapack_lp64 mkl_intel_lp64 mkl_sequential mkl_core mkl_blacs_intelmpi_lp64) 
target_link_libraries(pMat ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 

target_include_directories(pMat PUBLIC include ${MKLROOT}/include)


add_executable(parser_drive src/parser_drive.cpp)
target_link_libraries(parser_drive parser)





