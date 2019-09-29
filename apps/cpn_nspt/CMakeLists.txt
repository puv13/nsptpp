include_directories(../../src)

add_executable(cpn_nspt
  cpn_nspt_params.cpp cpn_nspt.cpp )


target_link_libraries(cpn_nspt
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
