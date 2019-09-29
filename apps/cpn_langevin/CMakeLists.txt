include_directories(../../src)

add_executable(cpn_langevin
  cpn_langevin_params.cpp cpn_langevin.cpp)


target_link_libraries(cpn_langevin
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
