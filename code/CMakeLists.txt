cmake_minimum_required(VERSION 3.0)

set(PROJECT graphene)

MESSAGE(STATUS "Using toolchain file: ${CMAKE_TOOLCHAIN_FILE}")

project(${PROJECT} C CXX)

#c++14 is required for some templates used in the code
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#find_package(Python 3.7 COMPONENTS Interpreter Development REQUIRED)
#if(Python_FOUND)
#    message(STATUS "Python Found: ${Python_EXECUTABLE}")
#    message(STATUS "Python Found: ${Python_INCLUDE_DIRS}")
#    message(STATUS "Python Found: ${Python_LIBRARIES}")
#    message(STATUS "Python Found: ${Python_LIBRARY_DIRS}")
#    include_directories(${Python_INCLUDE_DIRS})
#endif()

#FIND_PACKAGE(
#    Boost 1.79 REQUIRED
#    PATHS /home/ngolubev/Packages/libs/boost/1.80.0/include/
#)
#INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})

#Boost
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIB_DIR})

#Eigen
include_directories(${Eigen_PATH})

#mINI library
include_directories(${mINI_INCLUDE_DIR})

#cxxopts library
include_directories(${COPTS_INCLUDE_DIR})

#uBLAS odeint
include_directories("../libs/odeint/")

#SPLINTER
include_directories(${SPLINTER_INCLUDE_DIR})
link_directories(${SPLINTER_LIB_DIR})

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

find_package(OpenMP REQUIRED)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})

set(PROGS "main;analysis;data_writer;rspace")
#set(PROGS "KGMK")

foreach(PROG ${PROGS})
    set(BIN ${PROG}.exe)

    add_executable(${BIN} "${CMAKE_CURRENT_SOURCE_DIR}/${PROG}.cpp")

    target_link_libraries(${BIN} PUBLIC "stdc++fs")
    target_link_libraries(${BIN} PUBLIC ${SPLINTER_LIBS})

    target_link_libraries(${BIN} PUBLIC OpenMP::OpenMP_CXX)
ENDFOREACH()