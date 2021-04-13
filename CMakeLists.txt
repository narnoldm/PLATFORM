cmake_minimum_required(VERSION 3.16.2)
project(PDP VERSION 3.0)


find_package(MPI)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 17)
set(TECDIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_Party/tecio/teciosrc)
set(EIGENDIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_Party/eigen-3.3.7)


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
message("${CMAKE_CXX_COMPILER_ID} detected")
set(isINTEL 1)
message("MKLROOT $ENV{MKLROOT}")
set(MKLROOT $ENV{MKLROOT})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -DUSE_MKL")
endif()




add_subdirectory(src)
add_subdirectory(CASLAB_execs)
add_subdirectory(googletest)



#add_executable(parser_drive src/parser_drive.cpp)
#target_include_directories(parser_drive PUBLIC include)
#target_link_libraries(parser_drive parser meta)
#target_link_libraries(parser_drive debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)


#add_executable(serial_norm src/serial_norm_calc.cpp)
#target_include_directories(serial_norm PUBLIC include)
#target_link_libraries(serial_norm parser meta)
#target_link_libraries(serial_norm debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)








#enable_testing()


#add_executable(test_pMat tests/test_execs/pMat_test.cpp)
#target_include_directories(test_pMat PRIVATE include)
#target_link_libraries(test_pMat pMat)

#add_executable(test_meta tests/test_meta/test_bin.cpp)
#target_include_directories(test_meta PRIVATE include 3rd_Party/tecio/teciosrc)
#target_link_libraries(test_meta meta)

#add_executable(test_tec tests/test_meta/test_tecplot.cpp)
#target_include_directories(test_tec PRIVATE include 3rd_Party/tecio/teciosrc)
#target_link_libraries(test_tec meta)

#add_executable(test_tecSVD tests/test_meta/test_tecplotSVD.cpp)
#target_include_directories(test_tecSVD PRIVATE include 3rd_Party/tecio/teciosrc)
#target_link_libraries(test_tecSVD meta)

#add_executable(test_EvenRead tests/test_meta/test_EvenRead.cpp)
#target_include_directories(test_EvenRead PRIVATE include 3rd_Party/tecio/teciosrc)
#target_link_libraries(test_EvenRead meta)

#add_executable(test_EvenTec tests/test_meta/test_EvenTec.cpp)
#target_include_directories(test_EvenTec PRIVATE include 3rd_Party/tecio/teciosrc)
#target_link_libraries(test_EvenTec meta)


#add_executable(test_binSVD tests/test_meta/test_binSVD.cpp)
#target_include_directories(test_binSVD PRIVATE include 3rd_Party/tecio/teciosrc)
#target_link_libraries(test_binSVD meta)

#add_test(NAME test_pMat_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_pMat)
#add_test(NAME test_pMat_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_pMat)
#add_test(NAME test_metadata_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_meta)
#add_test(NAME test_metadata_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_meta)
#add_test(NAME test_EvenRead_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_EvenRead)
#add_test(NAME test_EvenRead_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_EvenRead)
#add_test(NAME test_EvenTec_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_EvenTec)
#add_test(NAME test_EvenTec_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_EvenTec)
#add_test(NAME test_EvenTec_parallel9 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 9 ${CMAKE_BINARY_DIR}/test_EvenTec)
#add_test(NAME test_tecplot_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_tec)
#add_test(NAME test_tecplot_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_tec)
#add_test(NAME test_tecplotSVD_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_tecSVD)
#add_test(NAME test_tecplotSVD_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_tecSVD)
#add_test(NAME test_binSVD_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_binSVD)
#add_test(NAME test_binSVD_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_binSVD)


#add_test(NAME test_SVD_bin_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/parser_drive)
#add_test(NAME test_SVD_bin_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/parser_drive)

#add_test(NAME test_cleanup WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND rm -r out even out1 load)








