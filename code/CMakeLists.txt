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

include_directories("/home/ngolubev/Packages/libs/boost/1.80.0/include/")
link_directories("/home/ngolubev/Packages/libs/boost/1.80.0/lib")

#FIND_PACKAGE(
#    Boost 1.79 REQUIRED
#    PATHS /home/ngolubev/Packages/libs/boost/1.80.0/include/
#)
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

#Eigen
set(Eigen_PATH "/home/ngolubev/Packages/libs/Eigen/eigen-3.4.0")
include_directories(${Eigen_PATH})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

#mINI library
include_directories("../libs/mINI/src/")

#uBLAS odeint
include_directories("../libs/odeint/")

#SPLINTER
include_directories("~/Packages/libs/splinter/include/")
link_directories("~/Packages/libs/splinter/build/")

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

#set(CMAKE_EXE_LINKER_FLAGS "-static")

find_package(OpenMP)

set(PROGS "main;analysis;data_writer;rspace")

foreach(PROG ${PROGS})
    set(BIN ${PROG}.exe)

    add_executable(${BIN} "${CMAKE_CURRENT_SOURCE_DIR}/${PROG}.cpp")

    target_link_libraries(${BIN} PUBLIC "stdc++fs")
    target_link_libraries(${BIN} PUBLIC splinter-static-3-0)

    if(OpenMP_CXX_FOUND)
        target_link_libraries(${BIN} PUBLIC OpenMP::OpenMP_CXX)
    endif()
ENDFOREACH()