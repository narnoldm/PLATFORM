



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(EX_1 test_param.cpp)
target_link_libraries(EX_1 meta param)

