find_package(MPI REQUIRED)

set(CMAKE_C_FLAGS          "${MPI_C_COMPILE_FLAGS} ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS        "${MPI_C_COMPILE_FLAGS} ${CMAKE_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${MPI_C_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")

# turn off inclusion of mpicxx.h
add_definitions(-DMPICH_SKIP_MPICXX)
add_definitions(-DOMPI_SKIP_MPICXX)

include_directories(${MPI_C_INCLUDE_PATH})

list(APPEND LIBRARIES ${MPI_C_LIBRARIES})
