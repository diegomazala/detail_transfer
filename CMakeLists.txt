cmake_minimum_required (VERSION 2.6)

project (DetailTransfer)

find_package( Eigen3 REQUIRED )
if( NOT EIGEN3_INCLUDE_DIR )
    set( EIGEN3_INCLUDE_DIR "$ENV{EIGEN3_INCLUDE_DIR}" )
endif()

find_package( OpenMesh REQUIRED )
find_package( OpenMesh REQUIRED )
if(NOT OPENMESH_FOUND) 
    message(Error "Could not find OpenMesh")
endif()


find_package(Threads)

add_subdirectory(src)


