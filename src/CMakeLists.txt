




set(SOURCES ${SOURCES})

include_directories(${PROJECT_SOURCE_DIR}/include)




add_library(pMat STATIC pMat.cpp processGrid.cpp)

add_dependencies(pMat Eigen3)

target_link_libraries(pMat PUBLIC ${MPI_C_LIBRARIES})
if(${isINTEL})
target_link_libraries(pMat PRIVATE ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 
target_include_directories(pMat PUBLIC ${MKLROOT}/include ${EIGEN_INSTALL_LOCATION}/include)
else()
target_link_libraries(pMat PRIVATE ${SCALAPACK_LIB})
endif()


add_library(meta STATIC metadata.cpp misc_func.cpp)
target_include_directories(meta PUBLIC ${TECDIR})
target_link_libraries(meta pMat ${TECDIR}/libtecio.a)
add_library(param STATIC param.cpp)
target_include_directories(param PRIVATE)