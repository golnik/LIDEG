set(CMAKE_FC_COMPILER "gfortran")
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")

#Boost
#set(Boost_INCLUDE_DIR "~/Packages/libs/boost/1.81.0/gcc8.3.0/include/")
#set(Boost_LIB_DIR "~/Packages/libs/boost/1.81.0/gcc8.3.0/lib/")
set(Boost_INCLUDE_DIR "/opt/ohpc/pub/apps/boost/1.79/include/")
set(Boost_LIB_DIR "/opt/ohpc/pub/apps/boost/1.79/lib/")

set(libs_PATH "/groups/ngolubev/Packages/libs/")

#Eigen
set(Eigen_PATH "${libs_PATH}/Eigen/eigen-3.4.0")

#SPLINTER
set(SPLINTER_INCLUDE_DIR "${libs_PATH}/splinter/include/")
set(SPLINTER_LIB_DIR "${libs_PATH}/splinter/build/")
set(SPLINTER_LIBS "splinter-static-3-0")

#mINI library
set(mINI_INCLUDE_DIR "${libs_PATH}/mINI/src/")

#cxxopts
set(COPTS_INCLUDE_DIR "${libs_PATH}/cxxopts/include/")
