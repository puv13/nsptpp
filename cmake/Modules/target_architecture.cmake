if(${TARGET_ARCHITECTURE} MATCHES "mic")

    if(NOT ${CMAKE_CXX_COMPILER_ID} MATCHES "Intel")
        message(FATAL_ERROR "Please set the C/C++ compilers to the Intel compiler (export CC=icc && export CXX=icpc), remove everything in your build directory, and rerun cmake")
    endif()

    if("$ENV{MPSS_SYSROOTS}" STREQUAL "")
        message(FATAL_ERROR "MPSS not found!")
    endif()

    if($ENV{MPSS_VERSION} MATCHES "[0-9]+\\.[0-9]+\\.[0-9]+")
        set(MPSS_VERSION $ENV{MPSS_VERSION})
    else()
        set(MPSS_VERSION $ENV{MPSS_VERSION}.0)
    endif()

    string(REPLACE "." "" MPSS_VERSION ${MPSS_VERSION})

    if(${MPSS_VERSION} LESS 370)
        ERROR("MPSS too old! Require at least Intel MPSS 3.7")
    endif()

    set(CMAKE_SYSTEM_NAME Linux)
    set(CMAKE_SYSTEM_PROCESSOR k1om)
    set(CMAKE_SYSTEM_VERSION 1)
    set(_CMAKE_TOOLCHAIN_PREFIX  x86_64-k1om-linux-)
    set(CMAKE_FIND_ROOT_PATH $ENV{MPSS_SYSROOTS}/k1om-mpss-linux)
    set(ARCH_FLAG "-mmic")

elseif(${TARGET_ARCHITECTURE} MATCHES "avx512")

    if(${CMAKE_CXX_COMPILER_ID}       MATCHES "Intel")
        set(ARCH_FLAG "-xMIC-AVX512")
    elseif(${CMAKE_CXX_COMPILER_ID}   MATCHES "GNU")
        set(ARCH_FLAG "-march=knl")
    elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
        set(ARCH_FLAG "-march=knl")
    endif()

else()

    set(ARCH_FLAG "-march=native")

endif()

set(CMAKE_C_FLAGS   "${ARCH_FLAG} ${CMAKE_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${ARCH_FLAG} ${CMAKE_CXX_FLAGS}")
