



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(POD POD_bin.cpp)
target_link_libraries(POD meta param)


