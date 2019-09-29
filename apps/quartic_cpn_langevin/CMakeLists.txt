include_directories(../../src)

add_executable(quartic_cpn_langevin
  quartic_cpn_langevin_params.cpp quartic_cpn_langevin.cpp )



target_link_libraries(quartic_cpn_langevin
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
