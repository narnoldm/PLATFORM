enable_testing()


add_executable(test_pMat test_execs/pMat_test.cpp)
target_include_directories(test_pMat PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(test_pMat pMat)

add_executable(test_meta test_meta/test_bin.cpp)
target_include_directories(test_meta PRIVATE ${CMAKE_SOURCE_DIR}/include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_meta meta)

add_executable(test_tec test_meta/test_tecplot.cpp)
target_include_directories(test_tec PRIVATE ${CMAKE_SOURCE_DIR}/include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_tec meta)

add_executable(test_tecSVD test_meta/test_tecplotSVD.cpp)
target_include_directories(test_tecSVD PRIVATE ${CMAKE_SOURCE_DIR}/include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_tecSVD meta)

add_executable(test_EvenRead test_meta/test_EvenRead.cpp)
target_include_directories(test_EvenRead PRIVATE ${CMAKE_SOURCE_DIR}/include 3rd_Party/tecio/teciosrc)
target_link_libraries(test_EvenRead meta)



add_test(NAME test_pMat_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/bin/test_pMat)
add_test(NAME test_pMat_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/test_pMat)
add_test(NAME test_metadata_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/bin/test_meta)
add_test(NAME test_metadata_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/test_meta)
add_test(NAME test_EvenRead_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/bin/test_EvenRead)
add_test(NAME test_EvenRead_parallel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/test_EvenRead)
add_test(NAME test_tecplot_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/bin/test_tec)
add_test(NAME test_tecplot_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/test_tec)
add_test(NAME test_tecplotSVD_serial WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 1 ${CMAKE_BINARY_DIR}/bin/test_tecSVD)
add_test(NAME test_tecplotSVD_paralllel WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND mpirun -np 4 ${CMAKE_BINARY_DIR}/bin/test_tecSVD)





add_test(NAME test_cleanup WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test_wd COMMAND rm -r out out1)
