include_directories(../../src)

add_executable(canonical_cpn_mc
  canonical_cpn_mc_params.cpp canonical_cpn_mc.cpp)


target_link_libraries(canonical_cpn_mc
  ${LIBRARIES}
  ${HDF5_LIBRARIES}
  ${BLAS_LIBRARIES}
  )
