include_directories(../../src)

add_executable(cpn_measurements
  cpn_measurements_params.cpp  cpn_measurements.cpp)



target_link_libraries(cpn_measurements 
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
