
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp.in ${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp @ONLY)

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

include_directories(${GMPXX_INCLUDE_DIR} ${GMP_INCLUDE_DIR})
include_directories(include)

add_library(objlib OBJECT ${SOURCES})

set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

add_library(Sapo SHARED $<TARGET_OBJECTS:objlib>)
add_library(SapoStatic STATIC $<TARGET_OBJECTS:objlib>)

set(PROJECT_LINK_LIBS ${PROJECT_LINK_LIBS} GMP::gmpxx GMP::gmp GLPK::glpk)

target_include_directories(Sapo PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_link_libraries(Sapo LINK_PUBLIC ${PROJECT_LINK_LIBS})
target_link_libraries(SapoStatic LINK_PUBLIC ${PROJECT_LINK_LIBS})

file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)
file(GLOB STL_HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/STL/*.h)

set_target_properties(SapoStatic PROPERTIES OUTPUT_NAME Sapo)

set(INCLUDE_INSTALL_DIR include )
set(LIB_INSTALL_DIR lib )

install(TARGETS Sapo DESTINATION ${LIB_INSTALL_DIR})
install(TARGETS SapoStatic DESTINATION ${LIB_INSTALL_DIR})
install(FILES ${HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR}/sapo)
install(FILES ${STL_HEADERS} DESTINATION ${INCLUDE_INSTALL_DIR}/sapo/STL)

include(CMakePackageConfigHelpers)
configure_package_config_file(SapoConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/SapoConfig.cmake
  INSTALL_DESTINATION ${LIB_INSTALL_DIR}/Sapo/cmake
  PATH_VARS INCLUDE_INSTALL_DIR)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/SapoConfig.cmake
  DESTINATION lib/Sapo/cmake )

if(CLANG_FORMAT_EXE)
   add_custom_target( format-lib-code COMMAND ${CLANG_FORMAT_EXE} -i
                             -style=file ${SOURCES} ${HEADERS} )
endif()

add_test(libSapo_compilation "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR}
	                   --target Sapo -j 4 -- )

add_subdirectory(tests)
