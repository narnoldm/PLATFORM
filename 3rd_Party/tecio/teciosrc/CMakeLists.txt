



if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING "CMake build type (Debug|Release)" FORCE)
endif ()


#FIND_PACKAGE(Boost)
#IF (Boost_FOUND)
#    INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIR})
#ELSE()
#  Set (Boost_INCLUDE_DIR /3rdpartylibs/boost/1.59.0/build/linux64-centos6.5/include)
#  message ("Warning:  Cannot find boost include directories.  Hardcoding to ${Boost_INCLUDE_DIR}")
#ENDIF()
include_directories(.)

IF (WIN32)
    set(BaseFlags "/EHsc /MP /wd\"4996\" /D MSWIN /D TP_PROJECT_USES_BOOST /D BOOST_ALL_NO_LIB /D MAKEARCHIVE /D NO_THIRD_PARTY_LIBS /D NO_ASSERTS /I\"${Boost_INCLUDE_DIR}\"")
ELSE ()
	set(BaseFlags "-DLINUX -DLINUX64 -DTP_PROJECT_USES_BOOST -DBOOST_ALL_NO_LIB -DMAKEARCHIVE -DNO_THIRD_PARTY_LIBS -DNO_ASSERTS -fmessage-length=0 -fPIC -isystem ${Boost_INCLUDE_DIR} -fvisibility=hidden -w")
ENDIF ()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BaseFlags}")
                    
FILE(GLOB mainFiles "*.cpp" "*.h")


add_library(tecio STATIC ${mainFiles})
add_dependencies(tecio Boost)
target_include_directories(tecio PUBLIC ${Boost_INCLUDE_DIR})
# Had to add this so building tecio care package within loci-chem works...
set_target_properties(tecio PROPERTIES LINKER_LANGUAGE CXX)
