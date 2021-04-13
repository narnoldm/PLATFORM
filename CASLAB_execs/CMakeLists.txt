



include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_executable(POD POD_tec.cpp)
target_include_directories(POD PUBLIC include)
target_link_libraries(POD meta)
#target_link_libraries(POD debug -Wl,--allow-multiple-definition,--undefined=malloc,--undefined=_ZdaPv /home/narnoldm/arm/forge/lib/64/libdmallocxx.a)


add_executable(projEr projEr.cpp)
target_include_directories(projEr PUBLIC include)
target_link_libraries(projEr param meta)


add_executable(projErP projErP.cpp)
target_include_directories(projErP PUBLIC include)
target_link_libraries(projErP param meta)

