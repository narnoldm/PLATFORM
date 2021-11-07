
include_directories(${PROJECT_SOURCE_DIR}/include)


# Timing executables

add_executable(PDP_time PDP_time.cpp)
target_link_libraries(PDP_time param meta)


