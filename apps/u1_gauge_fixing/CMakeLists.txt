include_directories(../../src)

add_executable(u1_gauge_fixing
  u1_gauge_fixing_params.cpp  u1_gauge_fixing.cpp)



target_link_libraries(u1_gauge_fixing
  ${LIBRARIES}
  ${BLAS_LIBRARIES}
  ${HDF5_LIBRARIES}
  )
