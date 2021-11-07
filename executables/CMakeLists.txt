



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(POD POD_tec.cpp)
target_link_libraries(POD meta param)

add_executable(projEr projEr.cpp)
target_link_libraries(projEr param meta)


add_executable(projErP projErP.cpp)
target_link_libraries(projErP param meta)


add_executable(QDEIM QDEIM_pre.cpp)
target_link_libraries(QDEIM param meta)


add_executable(calcSamplingError calcSamplingError.cpp)
target_link_libraries(calcSamplingError param meta)

