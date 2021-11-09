



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(EX_1 ex1_param.cpp)
target_link_libraries(EX_1 meta param)



add_executable(EX_2 ex2_pMat.cpp)
target_link_libraries(EX_2 meta param)


add_executable(EX_3 ex3_testbin.cpp)
target_link_libraries(EX_3 meta param)