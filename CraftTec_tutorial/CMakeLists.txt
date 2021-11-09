



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(EX_1 ex1_param.cpp)
target_link_libraries(EX_1 meta param)



add_executable(EX_2 ex2_pMat.cpp)
target_link_libraries(EX_2 meta param)


add_executable(EX_3 ex3_testbin.cpp)
target_link_libraries(EX_3 meta param)


add_executable(EX_4 ex4_testbinSVD.cpp)
target_link_libraries(EX_4 meta param)

add_executable(EX_5 ex5_testtecSVD.cpp)
target_link_libraries(EX_5 meta param)



add_executable(demo demo_binary.cpp)
target_link_libraries(demo meta param)