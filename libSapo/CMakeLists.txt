
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp.in ${CMAKE_CURRENT_SOURCE_DIR}/src/Version.cpp @ONLY)

file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

include_directories(include include/STL)

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

install(TARGETS Sapo DESTINATION lib)
install(TARGETS SapoStatic DESTINATION lib)
install(FILES ${HEADERS} DESTINATION include/sapo)
install(FILES ${STL_HEADERS} DESTINATION include/sapo/STL)

if(CLANG_FORMAT_EXE)
   add_custom_target( format-lib-code COMMAND ${CLANG_FORMAT_EXE} -i
                             -style=file ${SOURCES} ${HEADERS} )
endif()
