include_directories(../../src)

add_executable(principal_chiral_langevin
  principal_chiral_langevin_params.cpp principal_chiral_langevin.cpp)


target_link_libraries(principal_chiral_langevin
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
