






include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../include)

add_library(pMat STATIC pMat.cpp)
target_link_libraries(pMat ${MPI_C_LIBRARIES})
if(${isINTEL})
target_link_libraries(pMat ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 
target_include_directories(pMat PUBLIC include ${MKLROOT}/include)
else()
target_link_libraries(pMat ${SCALAPACK_LIB})
endif()

target_include_directories(pMat PUBLIC include ${EIGENDIR})

add_library(meta STATIC metadata.cpp misc_func.cpp)
target_include_directories(meta PUBLIC include ${TECDIR})
target_link_libraries(meta pMat ${TECDIR}/libtecio.a)
add_library(param STATIC param.cpp)
target_include_directories(param PRIVATE include)