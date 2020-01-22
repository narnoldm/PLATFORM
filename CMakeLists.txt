cmake_minimum_required(VERSION 3.16.2)
project(PDP VERSION 3.0)

message("MKLROOT $ENV{MKLROOT}")
set(MKLROOT $ENV{MKLROOT})
set(CMAKE_CXX_STANDARD 11)
set(TECDIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_Party/tecio/teciosrc)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
message("${CMAKE_CXX_COMPILER_ID} detected")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()



add_library(pMat STATIC src/pMat.cpp)
target_link_libraries(pMat ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 
target_include_directories(pMat PRIVATE include ${MKLROOT}/include)


add_library(parser STATIC src/parser.cpp)
target_include_directories(parser PRIVATE include)
target_link_libraries(parser pMat)


add_library(meta STATIC src/metadata.cpp)
target_include_directories(meta PUBLIC include 3rd_Party/tecio/teciosrc)
target_link_libraries(meta pMat ${TECDIR}/libtecio.a)

add_executable(parser_drive src/parser_drive.cpp)
target_include_directories(parser_drive PRIVATE include)
target_link_libraries(parser_drive parser pMat)



enable_testing()


add_executable(test_pMat tests/test_execs/pMat_test.cpp)
target_include_directories(test_pMat PRIVATE include)
target_link_libraries(test_pMat pMat)

add_executable(test_meta tests/test_meta/test_tecplot.cpp)
target_include_directories(test_meta PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_meta pMat meta)



add_test(NAME test_pMat_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_pMat)
add_test(NAME test_pMat_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_pMat)
add_test(NAME test_SVD_bin_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/parser_drive)
add_test(NAME test_SVD_bin_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/parser_drive)
add_test(NAME test_metadata_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_meta)
add_test(NAME test_metadata_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_meta)










