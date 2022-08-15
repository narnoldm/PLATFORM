



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(POD POD_tec.cpp)
target_link_libraries(POD meta param)

#add_executable(projEr projEr.cpp)
#target_link_libraries(projEr param meta)

#add_executable(projErP projErP.cpp)
#target_link_libraries(projErP param meta)

add_executable(calcSamplingError calcSamplingError.cpp)
target_link_libraries(calcSamplingError param meta)

add_executable(serial_norm serial_norm_calc.cpp)
target_link_libraries(serial_norm param meta)

# Timing executables
add_executable(PDP_time PDP_time.cpp)
target_link_libraries(PDP_time param meta)

# Subdirectories
add_subdirectory(hyper_reduction)
add_subdirectory(error_calcs)
