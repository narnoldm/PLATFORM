set(BINARY ${CMAKE_PROJECT_NAME}_tst)



set(SOURCES ${TEST_SOURCES})

file(GLOB_RECURSE TEST_SOURCES LIST_DIRECTORIES false *.h *.cpp)


add_executable(${BINARY} ${TEST_SOURCES})
add_test(NAME ${BINARY} COMMAND ${BINARY})

target_link_libraries(${BINARY} PUBLIC pMat meta param gtest)