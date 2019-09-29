include_directories(../../src)

add_executable(u1_measurements
  u1_measurements_params.cpp u1_measurements.cpp)


target_link_libraries(u1_measurements
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
