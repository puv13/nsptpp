if(${BUILD_CBLAS_VERSION} MATCHES "[Oo]penBLAS")
  message(STATUS "Looking for OpenBLAS lib")
  # Simple implementation using existing cmake module
  set(BLA_VENDOR OpenBLAS)
  set(BLAS_FIND_REQUIRED ON)
  include(FindBLAS)
  add_definitions(-DBLAS_AVAIL="TRUE")
elseif  (${BUILD_CBLAS_VERSION} MATCHES "mkl|MKL")
  message(STATUS "Looking for MKL lib")
  # Simple implementation using existing cmake module
  # set(BLA_VENDOR Intel)
  # set(BLAS_FIND_REQUIRED ON)
  # set(BLAS_DIR $ENV{MKL_LIBDIR})  
  # # include(FindBLAS)
  # find_path(BLAS_INCLUDE_DIRS mkl.h
  # 	    $ENV{MKL_BASE}/lib
  #   	)
  # add_definitions(-DBLAS_AVAIL="TRUE")
  # add_definitions(-DUSING_MKL="TRUE")
  # message( STATUS "BLAS found: (lib: ${BLAS_LIBRARIES})")
    message(STATUS "Looking for MKL lib")
  # Simple implementation using existing cmake module
  set(BLA_VENDOR All)
  set(BLAS_FIND_REQUIRED ON)
  include(FindBLAS)
  message( STATUS "BLAS found: (lib: ${BLAS_LIBRARIES})")
  add_definitions(-DUSE_MKL)
else ()
  message(STATUS "Looking for any BLAS lib")
  # Simple implementation using existing cmake module
  set(BLA_VENDOR All)
  set(BLAS_FIND_REQUIRED ON)
  include(FindBLAS)
  message( STATUS "BLAS found: (lib: ${BLAS_LIBRARIES})")	
endif()  
if (BLAS_FOUND)
   if  (${BUILD_CBLAS_VERSION} MATCHES "mkl|MKL")
       message(STATUS "MKL found")
       find_path(BLAS_INCLUDE_DIRS mkl.h
    	/usr/include
    	/usr/local/include
    	/usr/include/openblas
    	$ENV{MKL_INCDIR}
    	$ENV{BLAS_HOME}/include
    	)
   else ()
   	find_path(BLAS_INCLUDE_DIRS cblas.h
    	/usr/include
    	/usr/local/include
    	/usr/include/openblas
    	$ENV{BLAS_HOME}
    	$ENV{BLAS_HOME}/include
    	)
   endif()		
endif()

if(${BLAS_INCLUDE_DIRS} MATCHES "BLAS_INCLUDE_DIRS-NOTFOUND" )
  set(BLAS_FOUND FALSE)
  if  (${BUILD_CBLAS_VERSION} MATCHES "mkl|MKL")
      message(WARNING "Could not find mkl.h")
  else ()
       message(WARNING "Could not find cblas.h")
  endif()
else()
  message( STATUS "BLAS found: (lib: ${BLAS_LIBRARIES} include:  ${BLAS_INCLUDE_DIRS})" )
  add_definitions(-DBLAS_AVAIL="TRUE")
endif()
