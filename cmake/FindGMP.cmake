# Try to find the GMP libraries
# GMP_FOUND - system has GMP lib
# GMP_INCLUDE_DIR - the GMP include directory
# GMP_LIBRARY_DIR - Directory where the GMP libraries are located
# GMP_LIB - the GMP libraries
# GMPXX_INCLUDE_DIR - the GMPXX include directory
# GMPXX_LIBRARY_DIR - Directory where the GMPXX libraries are located
# GMPXX_LIB - the GMPXX libraries

find_path(GMP_INCLUDE_DIR
  NAMES gmp.h
  PATHS
  /usr/include/
  /usr/local/include/
  /opt/local/include/
  /opt/homebrew/include/
  DOC "The directory containing the GMP header files"
  )

find_path(GMPXX_INCLUDE_DIR
  NAMES gmpxx.h
  PATHS
  /usr/include/
  /usr/local/include/
  /opt/local/include/
  /opt/homebrew/include/
  DOC "The directory containing the GMP header files"
  )

find_path(GMP_LIBRARY_DIR
  NAMES "libgmp.so" "libgmp.dylib"
  PATHS
  /usr/lib/
  /usr/lib64/
  /usr/lib/x86_64-linux-gnu/
  /usr/local/lib/
  /opt/local/lib/
  /opt/homebrew/lib/
  DOC "Directory containing the GMP library"
  )

find_library(GMP_LIB NAMES gmp
  PATHS ENV GMP_LIBRARY_DIR
  DOC "Path to the GMP library"
  )

find_path(GMPXX_LIBRARY_DIR
  NAMES "libgmpxx.so" "libgmpxx.dylib"
  PATHS
  /usr/lib/
  /usr/lib64/
  /usr/lib/x86_64-linux-gnu/
  /usr/local/lib/
  /opt/local/lib/
  /opt/homebrew/lib/
  DOC "Directory containing the GMPXX library"
  )

find_library(GMPXX_LIB NAMES gmpxx
  PATHS ENV GMPXX_LIBRARY_DIR
  DOC "Path to the GMPXX library"
  )

# Set GMP_FIND_VERSION to 6.0.0 if no minimum version is specified
if(NOT GMP_FIND_VERSION)
  if(NOT GMP_FIND_VERSION_MAJOR)
    set(GMP_FIND_VERSION_MAJOR 6)
  endif()
  if(NOT GMP_FIND_VERSION_MINOR)
    set(GMP_FIND_VERSION_MINOR 0)
  endif()
  if(NOT GMP_FIND_VERSION_PATCH)
    set(GMP_FIND_VERSION_PATCH 0)
  endif()
  set(GMP_FIND_VERSION
    "${GMP_FIND_VERSION_MAJOR}.${GMP_FIND_VERSION_MINOR}.${GMP_FIND_VERSION_PATCH}")
endif()


if(GMP_INCLUDE_DIR AND GMP_LIB)

  # This program will fail to compile if GMP is too old.
  # We prefer to perform this "test" at compile-time to
  # avoid problems with e.g. try_run() during cross-compilation.
  file(WRITE ${PROJECT_BINARY_DIR}/gmp-version-check.cpp ""
  "#include <gmp.h>\n"
  "\n"
  "#define GMP_FIND_VERSION_MAJOR ${GMP_FIND_VERSION_MAJOR}\n"
  "#define GMP_FIND_VERSION_MINOR ${GMP_FIND_VERSION_MINOR}\n"
  "#define GMP_FIND_VERSION_PATCH ${GMP_FIND_VERSION_PATCH}\n"
  "\n"
  "#if __GNU_MP_VERSION < GMP_FIND_VERSION_MAJOR\n"
  "#error insufficient GMP major version\n"
  "#elif __GNU_MP_VERSION == GMP_FIND_VERSION_MAJOR\n"
  "#if __GNU_MP_VERSION_MINOR < GMP_FIND_VERSION_MINOR\n"
  "#error insufficient GMP minor version\n"
  "#elif __GNU_MP_VERSION_MINOR == GMP_FIND_VERSION_MINOR\n"
  "#if __GNU_MP_VERSION_PATCH < GMP_FIND_VERSION_PATCH\n"
  "#error insufficient GMP patch version\n"
  "#endif\n"
  "#endif\n"
  "#endif\n"
  "#pragma GCC diagnostic ignored \"-Wunused-parameter\"\n"
  "\n"
  "int main(int argc, char** argv) { return 0; }\n")

  # Try to compile the test program above with the appropriate version
  # strings substituted in.
  try_compile(GMP_VERSION_OK
          "${PROJECT_BINARY_DIR}"
          "${PROJECT_BINARY_DIR}/gmp-version-check.cpp"
          CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${GMP_INCLUDE_DIR}")
endif()

if(NOT GMP_VERSION_OK)
  message(STATUS "No sufficient GMP version detected")
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(GMP
  DEFAULT_MSG
  GMPXX_INCLUDE_DIR GMPXX_LIB GMP_INCLUDE_DIR GMP_LIB GMP_VERSION_OK
)

# Set targets
if(GMP_FOUND)
  # C library
  if(NOT TARGET GMP::gmp)
    add_library(GMP::gmp UNKNOWN IMPORTED)
    set_target_properties(GMP::gmp PROPERTIES
      IMPORTED_LOCATION ${GMP_LIB}
      INTERFACE_INCLUDE_DIRECTORIES ${GMP_INCLUDE_DIR}
    )
  endif()

  # C++ library, which requires a link to the C library
  if(NOT TARGET GMP::gmpxx)
    add_library(GMP::gmpxx UNKNOWN IMPORTED)
    set_target_properties(GMP::gmpxx PROPERTIES
      IMPORTED_LOCATION ${GMPXX_LIB}
      INTERFACE_INCLUDE_DIRECTORIES ${GMPXX_INCLUDE_DIR}
      INTERFACE_LINK_LIBRARIES GMP::gmp
    )
  endif()
endif()

