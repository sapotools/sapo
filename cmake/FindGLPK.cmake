# This module defines the following uncached variables:
#  GLPK_FOUND, if false, do not try to use GLPK.
#  GLPK_THREADS, if false GLPK does not support thread local memory.
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

if(GLPK_INCLUDE_DIR AND GLPK_LIB)

  # Set GLPK_FIND_VERSION to 6.0.0 if no minimum version is specified
  if(NOT GLPK_FIND_VERSION)
  if(NOT GLPK_FIND_VERSION_MAJOR)
    set(GLPK_FIND_VERSION_MAJOR 5)
  endif()
  if(NOT GLPK_FIND_VERSION_MINOR)
    set(GLPK_FIND_VERSION_MINOR 0)
  endif()
  set(GLPK_FIND_VERSION
    "${GLPK_FIND_VERSION_MAJOR}.${GLPK_FIND_VERSION_MINOR}")
  endif()

  # This program will fail to compile if GLPK is too old and
  # fail to run if the GLPK library does not support thread 
  # local memory
  file(WRITE ${PROJECT_BINARY_DIR}/glpk-check.cpp ""
  "#include <glpk.h>\n"
  "#include <stdlib.h>\n"
  "\n"
  "#define GLPK_FIND_VERSION_MAJOR ${GLPK_FIND_VERSION_MAJOR}\n"
  "#define GLPK_FIND_VERSION_MINOR ${GLPK_FIND_VERSION_MINOR}\n"
  "\n"
  "#if GLP_MAJOR_VERSION < GLPK_FIND_VERSION_MAJOR\n"
  "#error insufficient GLPK major version\n"
  "#elif GLP_MINOR_VERSION == GLPK_FIND_VERSION_MINOR\n"
  "#if GLP_VERSION_MINOR < GLPK_FIND_VERSION_MINOR\n"
  "#error insufficient GLPK minor version\n"
  "#endif\n"
  "#endif\n"
  "#pragma GCC diagnostic ignored \"-Wunused-parameter\"\n"
  "\n"
  "int main(int argc, char** argv) { \n"
  " if(!glp_config(\"TLS\")) { exit(EXIT_FAILURE); } \n"
  " return 0;\n"
  "}\n")

  try_run(GLPK_THREADS GLPK_VERSION_OK
          "${PROJECT_BINARY_DIR}"
          "${PROJECT_BINARY_DIR}/glpk-check.cpp"
          CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${GLPK_INCLUDE_DIR}"
          LINK_LIBRARIES ${GLPK_LIB})
  
  if(GLPK_THREADS EQUAL 0)
    set(GLPK_THREADS TRUE)
  else()
    set(GLPK_THREADS FALSE)
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GLPK
  DEFAULT_MSG
  GLPK_INCLUDE_DIR GLPK_LIBRARY_DIR GLPK_LIB GLPK_VERSION_OK
)

if(GLPK_FOUND AND GLPK_VERSION_OK)
  # C library
  if(NOT TARGET GLPK::glpk)
    add_library(GLPK::glpk UNKNOWN IMPORTED)
    set_target_properties(GLPK::glpk PROPERTIES
      IMPORTED_LOCATION ${GLPK_LIB}
      INTERFACE_INCLUDE_DIRECTORIES ${GLPK_INCLUDE_DIR}
    )
  endif()
endif()