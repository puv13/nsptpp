include_directories(../../src)

add_executable(quartic_cpn_nspt
  quartic_cpn_nspt_params.cpp quartic_cpn_nspt.cpp )


target_link_libraries(quartic_cpn_nspt
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
