



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(POD POD_tec.cpp)
target_link_libraries(POD meta param)

add_executable(serial_norm_calc serial_norm_calc.cpp)
target_link_libraries(serial_norm_calc meta param)