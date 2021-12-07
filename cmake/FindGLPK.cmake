# This module defines the following uncached variables:
#  GLPK_FOUND, if false, do not try to use GLPK.
#  GLPK_INCLUDE_DIR, where to find glpk.h.
#  GLPK_LIB, the libraries to link against to use the GLPK library
#  GLPK_LIBRARY_DIR, the directory where the GLPK library is found.

find_path(GLPK_INCLUDE_DIR
  glpk.h
  PATHS
  /usr/include/
  /usr/local/include/
  /opt/local/include/
  /opt/homebrew/include/
  DOC "The directory containing the GLPK header files"
)

  find_library(GLPK_LIBRARY_DIR
    NAMES "libglpk.so" "libglpk.dylib" "libglpk.a"
    PATHS
    /usr/lib/
    /usr/lib64/
    /usr/lib/x86_64-linux-gnu/
    /usr/local/lib/
    /opt/local/lib/
    /opt/homebrew/lib/
    DOC "Directory containing the GLPK library"
  )

  find_library(GLPK_LIB NAMES glpk
   PATHS ENV GLPK_LIBRARY_DIR
   DOC "Path to the GLPK library"
  )
	    
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK
  DEFAULT_MSG
  GLPK_INCLUDE_DIR GLPK_LIBRARY_DIR GLPK_LIB
)

if(GLPK_FOUND)
  # C library
  if(NOT TARGET GLPK::glpk)
    add_library(GLPK::glpk UNKNOWN IMPORTED)
    set_target_properties(GLPK::glpk PROPERTIES
      IMPORTED_LOCATION ${GLPK_LIB}
      INTERFACE_INCLUDE_DIRECTORIES ${GLPK_INCLUDE_DIR}
    )
  endif()
endif()