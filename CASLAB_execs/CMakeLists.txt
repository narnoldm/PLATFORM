



include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(POD POD_tec.cpp)
target_include_directories(POD PUBLIC)
target_link_libraries(POD meta param)
#target_link_libraries(POD debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)


add_executable(projEr projEr.cpp)
target_include_directories(projEr PUBLIC)
target_link_libraries(projEr param meta)


add_executable(projErP projErP.cpp)
target_include_directories(projErP PUBLIC)
target_link_libraries(projErP param meta)


add_executable(QDEIM QDEIM_pre.cpp)
target_include_directories(QDEIM PUBLIC include)
target_link_libraries(QDEIM param meta)
#target_link_libraries(QDEIM debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)



add_executable(calcSamplingError calcSamplingError.cpp)
target_include_directories(calcSamplingError PUBLIC include)
target_link_libraries(calcSamplingError param meta)

# Timing executables

add_executable(PDP_time PDP_time.cpp)
target_include_directories(PDP_time PUBLIC include)
target_link_libraries(PDP_time param meta)