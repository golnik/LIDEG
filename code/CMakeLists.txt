cmake_minimum_required(VERSION 3.0)

set(PROJECT graphene)

set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

project(${PROJECT} C CXX)

#c++14 is required for some templates used in the code
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Python 3.7 COMPONENTS Interpreter Development REQUIRED)
if(Python_FOUND)
    message(STATUS "Python Found: ${Python_EXECUTABLE}")
    message(STATUS "Python Found: ${Python_INCLUDE_DIRS}")
    message(STATUS "Python Found: ${Python_LIBRARIES}")
    message(STATUS "Python Found: ${Python_LIBRARY_DIRS}")
    include_directories(${Python_INCLUDE_DIRS})
endif()

include_directories("/home/ngolubev/Packages/libs/boost/1.79.0/boost_1_79_0")

link_directories("/home/ngolubev/Packages/libs/boost/1.79.0/boost_1_79_0/stage/lib")

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#mINI library
include_directories("../libs/mINI/src/")

#uBLAS odeint
include_directories("../libs/odeint/")

#add_library(
#    ${PROJECT} SHARED
#    "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
#    "${CMAKE_CURRENT_SOURCE_DIR}/python_binding.cpp"
#)

#target_link_libraries(${PROJECT}
#    boost_python37
#    boost_numpy37
#    ${Python_LIBRARIES}
#)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
#set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

add_executable(
    main.exe
    "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"  
)

add_executable(
    analysis.exe
    "${CMAKE_CURRENT_SOURCE_DIR}/analysis.cpp"  
)

add_executable(
    rspace.exe
    "${CMAKE_CURRENT_SOURCE_DIR}/rspace.cpp"  
)

#add_executable(
#    debug.exe
#    "${CMAKE_CURRENT_SOURCE_DIR}/utils/gauss_legendre.c"    
#    "${CMAKE_CURRENT_SOURCE_DIR}/debug.cpp"
#)

target_link_libraries(main.exe PRIVATE "stdc++fs")
target_link_libraries(analysis.exe PRIVATE "stdc++fs")
target_link_libraries(rspace.exe PRIVATE "stdc++fs")
#target_link_libraries(debug.exe PRIVATE "stdc++fs")

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(main.exe PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(analysis.exe PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(rspace.exe PUBLIC OpenMP::OpenMP_CXX)
endif()