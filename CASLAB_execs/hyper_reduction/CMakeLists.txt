
add_library(sampling STATIC src/sampling.cpp)
target_include_directories(sampling PUBLIC ${PROJECT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR}/CASLAB_execs/hyper_reduction/include)
target_link_libraries(sampling pMat param)

add_executable(QDEIM QDEIM_pre.cpp)
target_link_libraries(QDEIM param meta sampling)
