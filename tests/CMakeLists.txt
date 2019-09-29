# No longer needed. We build gtest ourselves 
# find_package(GTest REQUIRED)

find_package(Threads REQUIRED)

include_directories(${GTEST_INCLUDE_DIRS})

include_directories(../src)

add_executable(alltests
    EXCLUDE_FROM_ALL
    test_expansion.cpp
    test_su3.cpp
    test_lattice.cpp	
    test_io.cpp
    test_cp.cpp
    test_qcd.cpp
    test_u1.cpp
    test_canonical_cpn.cpp
    test_quartic_cpn.cpp
    test_compactQED.cpp
    test_nummat.cpp
    test_sun.cpp
    test_randgen.cpp
    test_stat.cpp
    test_pcm.cpp
    )

target_link_libraries(alltests
    ${LIBRARIES}
    # No longer defined 
    #${GTEST_BOTH_LIBRARIES}
    ${CMAKE_THREAD_LIBS_INIT}
    ${HDF5_LIBRARIES}
    ${BLAS_LIBRARIES}
    # Tell cmake about gtest dependency 
    gtest
    gtest_main
    )

add_test(NAME test COMMAND alltests)
install( PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/alltests DESTINATION tests/nsptpp/ )
