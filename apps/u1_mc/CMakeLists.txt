include_directories(../../src)

add_executable(u1_mc
  u1_mc_params.cpp u1_mc.cpp)


target_link_libraries(u1_mc
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
