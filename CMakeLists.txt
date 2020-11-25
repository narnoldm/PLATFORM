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



add_library(pMat STATIC src/pMat.cpp)
target_link_libraries(pMat ${MPI_C_LIBRARIES})
if(${isINTEL})
target_link_libraries(pMat ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 
target_include_directories(pMat PUBLIC include ${MKLROOT}/include)
else()
target_link_libraries(pMat -L/home/narnoldm/Dropbox/Michigan/Lab/0_2020/1Jan/ScaLAPACK/build/lib/libscalapack.a)
endif()
target_include_directories(pMat PUBLIC include ${EIGENDIR})
#target_link_libraries(pMat debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)

add_library(meta STATIC src/metadata.cpp src/misc_func.cpp)
target_include_directories(meta PUBLIC include ${TECDIR})
target_link_libraries(meta pMat ${TECDIR}/libtecio.a)
#target_link_libraries(meta debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)



add_library(parser STATIC src/dataID.cpp src/executioner.cpp src/inputInfo.cpp src/inputReader.cpp src/misc_func.cpp src/operation.cpp src/sysInfo.cpp)
target_include_directories(parser PRIVATE include)
target_link_libraries(parser meta)
#target_link_libraries(parser debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)

add_library(param STATIC src/param.cpp)
target_include_directories(param PRIVATE include)


add_executable(parser_drive src/parser_drive.cpp)
target_include_directories(parser_drive PUBLIC include)
target_link_libraries(parser_drive parser meta)
#target_link_libraries(parser_drive debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)


add_executable(serial_norm CASLAB_execs/serial_norm_calc.cpp)
target_include_directories(serial_norm PUBLIC include)
target_link_libraries(serial_norm param meta)
#target_link_libraries(serial_norm debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)


add_executable(POD CASLAB_execs/POD_tec.cpp)
target_include_directories(POD PUBLIC include)
target_link_libraries(POD parser param meta)
#target_link_libraries(POD debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)

add_executable(QDEIM CASLAB_execs/QDEIM_pre.cpp)
target_include_directories(QDEIM PUBLIC include)
target_link_libraries(QDEIM parser param meta)
#target_link_libraries(QDEIM debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)


add_executable(projEr CASLAB_execs/projEr.cpp)
target_include_directories(projEr PUBLIC include)
target_link_libraries(projEr parser param meta)


add_executable(projErP CASLAB_execs/projErP.cpp)
target_include_directories(projErP PUBLIC include)
target_link_libraries(projErP parser param meta)


add_executable(calcSamplingError CASLAB_execs/calcSamplingError.cpp)
target_include_directories(calcSamplingError PUBLIC include)
target_link_libraries(calcSamplingError parser param meta)

# Timing executables

add_executable(PDP_time tests/timing/PDP_time.cpp)
target_include_directories(PDP_time PUBLIC include)
target_link_libraries(PDP_time parser param meta)



enable_testing()

add_executable(test_param tests/test_execs/test_param.cpp)
target_include_directories(test_param PRIVATE include)
target_link_libraries(test_param param ${MPI_C_LIBRARIES})

add_executable(test_pMat tests/test_execs/pMat_test.cpp)
target_include_directories(test_pMat PRIVATE include)
target_link_libraries(test_pMat pMat)

add_executable(test_pinv tests/test_execs/test_pinv.cpp)
target_include_directories(test_pinv PRIVATE include)
target_link_libraries(test_pinv pMat)


add_executable(test_meta tests/test_meta/test_bin.cpp)
target_include_directories(test_meta PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_meta meta)

add_executable(test_tec tests/test_meta/test_tecplot.cpp)
target_include_directories(test_tec PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_tec meta)

add_executable(test_tecSVD tests/test_meta/test_tecplotSVD.cpp)
target_include_directories(test_tecSVD PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_tecSVD meta)

add_executable(test_EvenRead tests/test_meta/test_EvenRead.cpp)
target_include_directories(test_EvenRead PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_EvenRead meta)

add_executable(test_EvenTec tests/test_meta/test_EvenTec.cpp)
target_include_directories(test_EvenTec PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_EvenTec meta)


add_executable(test_binSVD tests/test_meta/test_binSVD.cpp)
target_include_directories(test_binSVD PRIVATE include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_binSVD meta)


add_test(NAME test_param_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_param)
add_test(NAME test_param_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_param)
add_test(NAME test_pMat_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_pMat)
add_test(NAME test_pMat_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_pMat)
add_test(NAME test_pinv_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_pinv)
add_test(NAME test_pinv_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_pinv)
add_test(NAME test_metadata_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_meta)
add_test(NAME test_metadata_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_meta)
add_test(NAME test_EvenRead_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_EvenRead)
add_test(NAME test_EvenRead_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_EvenRead)
add_test(NAME test_EvenTec_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_EvenTec)
add_test(NAME test_EvenTec_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_EvenTec)
add_test(NAME test_EvenTec_parallel9 WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 9 ${CMAKE_BINARY_DIR}/test_EvenTec)
add_test(NAME test_tecplot_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_tec)
add_test(NAME test_tecplot_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_tec)
add_test(NAME test_tecplotSVD_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_tecSVD)
add_test(NAME test_tecplotSVD_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_tecSVD)
add_test(NAME test_binSVD_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/test_binSVD)
add_test(NAME test_binSVD_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/test_binSVD)


add_test(NAME test_SVD_bin_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/parser_drive)
add_test(NAME test_SVD_bin_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test1_SVD_bin COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/parser_drive)



add_test(NAME test_cleanup WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_wd COMMAND rm -r out even load out1 )








