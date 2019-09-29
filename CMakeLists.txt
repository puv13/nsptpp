# The MIT License (MIT)

# Copyright (c) 2019 Jakob Simeth and Matthias Puhr
# Copyright (c) 2015 Crascit

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.



cmake_minimum_required(VERSION 3.1)

set(PROJECT "nsptpp")
project(${PROJECT} VERSION "0.1.0")

enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set(TARGET_ARCHITECTURE  "x86" CACHE STRING "build for target architecture (avx512|mic|x86)")
set(BUILD_OPENMP_VERSION ON    CACHE BOOL   "build version with OpenMP")
set(BUILD_MPI_VERSION    OFF   CACHE BOOL   "build version with MPI")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#link_directories($ENV{MKL_LIBDIR}) 

#Set Linker flags
#set(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")


# ######################################################################
# GTEST
# ######################################################################

# This is a modification of code from
# https://crascit.com/2015/07/25/cmake-gtest/
# See also https://github.com/Crascit/DownloadProject
# for a more sophisticated version. (MIT license)


# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
execute_process(COMMAND "${CMAKE_COMMAND}" --build .
  WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/googletest-download" )
 
# Prevent GoogleTest from overriding our compiler/linker options
# when building with Visual Studio
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This adds
# the following targets: gtest, gtest_main, gmock
# and gmock_main
add_subdirectory("${CMAKE_BINARY_DIR}/googletest-src"
                 "${CMAKE_BINARY_DIR}/googletest-build")
 
########################################################################

include(target_architecture)
include(language_standard)
include(optimization_level)
include(warnings)

if(${BUILD_OPENMP_VERSION})
    include(openmp)
endif()

if(${BUILD_MPI_VERSION})
    include(mpi)
endif()

# intentionally without ${} !
#if(BUILD_CBLAS_VERSION)
include(cblas)
if (BLAS_FOUND)
  message(STATUS "Blas include: ${BLAS_INCLUDE_DIRS}.")
  include_directories(${BLAS_INCLUDE_DIRS})
else()
  message(WARNING "Building without BLAS support.")
endif()

  

execute_process(COMMAND git rev-parse --verify HEAD
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
    )
string(TIMESTAMP BUILD_TIMESTAMP)

configure_file( ${CMAKE_SOURCE_DIR}/src/nsptpp_config.h.in
  "${CMAKE_CURRENT_BINARY_DIR}/nsptpp_config.h" @ONLY)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

#-- handle  highfive
set(HIGHFIVE "otherlibs/HighFive")
include_directories(${HIGHFIVE}/include)
install(DIRECTORY ${HIGHFIVE}/include/ DESTINATION include/nsptpp)

find_package(HDF5 REQUIRED COMPONENTS C HL)
include_directories(${HDF5_INCLUDE_DIRS})



#add_subdirectory(src)

add_subdirectory(tests)
add_subdirectory(apps/canonical_cpn_mc)
add_subdirectory(apps/cpn_measurements)
add_subdirectory(apps/cpn_gauge_fixing)
add_subdirectory(apps/u1_mc)
add_subdirectory(apps/u1_gauge_fixing)
add_subdirectory(apps/u1_measurements)
add_subdirectory(apps/cpn_langevin)
add_subdirectory(apps/cpn_nspt)
add_subdirectory(apps/quartic_cpn_nspt)
add_subdirectory(apps/quartic_cpn_langevin)
add_subdirectory(apps/principal_chiral_langevin)
add_subdirectory(apps/principal_chiral_nspt)

