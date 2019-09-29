include_directories(../../src)

add_executable(principal_chiral_nspt
  principal_chiral_nspt_params.cpp principal_chiral_nspt.cpp)


target_link_libraries(principal_chiral_nspt
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
