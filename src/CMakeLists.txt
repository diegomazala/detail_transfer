cmake_minimum_required (VERSION 2.6)


set (TargetName detail_transfer)
project (${TargetName})
add_executable(${TargetName} 
	detail_displace.cpp 
	tinyply.cpp 
	tinyply.h)
target_compile_features(${TargetName} PRIVATE cxx_std_17)
target_include_directories(${TargetName} PRIVATE ${OPENMESH_INCLUDE_DIRS})
target_link_libraries(${TargetName} PRIVATE ${OPENMESH_LIBRARIES} stdc++fs Threads::Threads)

set (TargetName smooth)
project (${TargetName})
add_executable(${TargetName} 
	smooth.cpp 
	smooth.h)
add_definitions(-D_USE_MATH_DEFINES)
target_compile_features(${TargetName} PRIVATE cxx_std_11)
target_include_directories(${TargetName} PRIVATE ${OPENMESH_INCLUDE_DIRS})
target_link_libraries(${TargetName} PRIVATE ${OPENMESH_LIBRARIES} stdc++fs)

