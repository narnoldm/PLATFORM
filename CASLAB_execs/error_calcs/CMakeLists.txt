
add_library(error STATIC src/error_funcs.cpp)
target_include_directories(error PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/CASLAB_execs/error_calcs/include)
target_link_libraries(error pMat meta)

add_executable(field_error field_error.cpp)
target_link_libraries(field_error param meta error)

add_executable(sampling_error sampling_error.cpp)
target_link_libraries(sampling_error param meta error)