include_directories(../../src)

add_executable(cpn_gauge_fixing
  cpn_gauge_fixing_params.cpp  cpn_gauge_fixing.cpp)



target_link_libraries(cpn_gauge_fixing
  ${LIBRARIES}
  ${BLAS_LIBRARIES}
  ${HDF5_LIBRARIES}
  )
