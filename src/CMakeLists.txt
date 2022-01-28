




set(SOURCES ${SOURCES})

include_directories(${PROJECT_SOURCE_DIR}/include)




add_library(pMat STATIC pMat.cpp processGrid.cpp)

add_dependencies(pMat Eigen3)


if(${isINTEL})
target_link_libraries(pMat PRIVATE ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a -Wl,--end-group -lpthread -lm -ldl) 
target_include_directories(pMat PUBLIC ${MKLROOT}/include ${EIGEN_INSTALL_LOCATION}/include)
else()
add_dependencies(pMat Scalapack)
message("${SCALAPACK_LIB_DIR}")
target_include_directories(pMat PUBLIC)
target_link_libraries(pMat PUBLIC scalapack)
endif()
target_link_libraries(pMat PUBLIC ${MPI_C_LIBRARIES})
target_include_directories(pMat PUBLIC ${Eigen_INCLUDE_DIR})

add_library(meta STATIC metadata.cpp misc_func.cpp misc_func_mpi.cpp)
target_include_directories(meta PUBLIC ${TECDIR})
target_link_libraries(meta pMat tecio)
add_library(param STATIC param.cpp)
target_include_directories(param PRIVATE)