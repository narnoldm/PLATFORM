cmake_minimum_required(VERSION 3.16.2)
project(PDP VERSION 3.0)

message("MKLROOT $ENV{MKLROOT}")
set(MKLROOT $ENV{MKLROOT})
set(CMAKE_CXX_STANDARD 17)
set(TECDIR ${CMAKE_CURRENT_SOURCE_DIR}/3rd_Party/tecio/teciosrc)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
message("${CMAKE_CXX_COMPILER_ID} detected")
set(isINTEL 1)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
endif()



add_library(pMat STATIC src/pMat.cpp)
target_link_libraries(pMat ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 
endif()
target_include_directories(pMat PUBLIC include ${MKLROOT}/include)
target_link_libraries(pMat debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/allinea/forge/lib/64/libdmallocxx.a)

add_library(meta STATIC src/metadata.cpp)
target_include_directories(meta PUBLIC include 3rd_Party/tecio/teciosrc)
target_link_libraries(meta pMat ${TECDIR}/libtecio.a)
target_link_libraries(meta debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/allinea/forge/lib/64/libdmallocxx.a)



add_library(parser STATIC src/dataID.cpp src/executioner.cpp src/inputInfo.cpp src/inputReader.cpp src/misc_func.cpp src/operation.cpp src/sysInfo.cpp)
target_include_directories(parser PRIVATE include)
target_link_libraries(parser meta)
target_link_libraries(parser debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/allinea/forge/lib/64/libdmallocxx.a)


add_executable(parser_drive src/parser_drive.cpp)
target_include_directories(parser_drive PUBLIC include)
target_link_libraries(parser_drive parser meta)
target_link_libraries(parser_drive debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/allinea/forge/lib/64/libdmallocxx.a)



enable_testing()


add_executable(test_pMat tests/test_execs/pMat_test.cpp)
target_include_directories(test_pMat PRIVATE include)
target_link_libraries(test_pMat pMat)

add_executable(test_meta tests/test_meta/test_bin.cpp)
target_include_directories(test_meta PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_meta meta)

add_executable(test_tec tests/test_meta/test_tecplot.cpp)
target_include_directories(test_tec PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_tec meta)

add_executable(test_tecSVD tests/test_meta/test_tecplotSVD.cpp)
target_include_directories(test_tecSVD PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_tecSVD meta)

add_test(NAME test_pMat_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_pMat)
add_test(NAME test_pMat_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_pMat)
add_test(NAME test_metadata_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_meta)
add_test(NAME test_metadata_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_meta)
add_test(NAME test_tecplot_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_tec)
add_test(NAME test_tecplot_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_tec)
add_test(NAME test_tecplotSVD_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_tecSVD)
add_test(NAME test_tecplotSVD_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_tecSVD)
add_test(NAME test_SVD_bin_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/parser_drive)
add_test(NAME test_SVD_bin_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/parser_drive)



add_test(NAME test_cleanup WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND rm -r out)








