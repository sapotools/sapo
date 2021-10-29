cmake_minimum_required(VERSION 3.7)

project (sapo CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#For shared libraries:
set ( PROJECT_LINK_LIBS ginac glpk cln )

include_directories(include include/STL parser/include)

find_package(FLEX 2.6.3 REQUIRED)
find_package(BISON REQUIRED)

flex_target(input_scanner parser/scanner.l ${CMAKE_SOURCE_DIR}/parser/src/scanner.cpp)
bison_target(input_parser parser/parser.yy ${CMAKE_SOURCE_DIR}/parser/src/parser.cpp
							DEFINES_FILE ${CMAKE_SOURCE_DIR}/parser/include/parser.h)

# This is to support brew
include_directories(/usr/local/Cellar/ginac/1.8.1/include) 
include_directories(/usr/local/Cellar/cln/1.3.4/include)
include_directories(/usr/local/Cellar/glpk/5.0/include)
include_directories(⁨/usr/local/include) 

link_directories( /usr/local/Cellar/ginac/1.8.1/lib )
link_directories( /usr/local/Cellar/cln/1.3.4/lib )
link_directories( /usr/local/Cellar/glpk/5.0/lib )
link_directories( /usr/local/lib )
# end brew support


file(GLOB_RECURSE SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp ${CMAKE_SOURCE_DIR}/parser/src/*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_executable(sapo ${BISON_input_parser_OUTPUT_SOURCE} ${FLEX_input_scanner_OUTPUTS} ${SOURCES})

set(CMAKE_CXX_FLAGS "-O2 -Wall -Werror")

target_compile_features(sapo PRIVATE cxx_range_for)
target_link_libraries(sapo ${PROJECT_LINK_LIBS} )

find_program( CLANG_FORMAT_EXE NAMES "clang-format")
if(NOT CLANG_FORMAT_EXE)
   message(STATUS "clang-format not found.")
else()
   message(STATUS "clang-format found: ${CLANG_FORMAT_EXE}")
endif()

if(CLANG_FORMAT_EXE)
   file(GLOB_RECURSE HEADERS ${CMAKE_SOURCE_DIR}/include/*.h ${CMAKE_SOURCE_DIR}/parser/include/*.h)

   add_custom_target( format-code COMMAND ${CLANG_FORMAT_EXE} -i -style=file ${SOURCES} ${HEADERS} )
endif()

execute_process(COMMAND sh ${CMAKE_SOURCE_DIR}/scripts/install-hooks.bash )

enable_testing()

set(sapo_tests Ebola Influenza Influenza-no-splits Influenza-splits LotkaVolterra
               Phosphorelay Quadcopter Rossler SIRp SIR VanDerPol )

add_test(sapo_compilation "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR} --target sapo -j 4 -- )
foreach(test_name ${sapo_tests})
   add_test(NAME ${test_name}
	    COMMAND ${CMAKE_COMMAND} -DTEST_NAME=${test_name} -DTEST_DIR=${CMAKE_SOURCE_DIR}/examples
			 -P ${CMAKE_SOURCE_DIR}/runtests.cmake)
   set_tests_properties(${test_name} PROPERTIES FIXTURES_SETUP sapo_compilation)
endforeach()
