cmake_minimum_required(VERSION 2.8.9)
project (sapo CXX)

#For shared libraries:
set ( PROJECT_LINK_LIBS ginac glpk)
link_directories( /usr/local/lib )

find_package(FLEX)
find_package(BISON)

flex_target(input_scanner parser/scanner.l parser/src/scanner.cpp)
bison_target(input_parser parser/parser.yy parser/src/parser.cpp
							DEFINES_FILE parser/include/parser.h)

include_directories(⁨/usr/local/include /usr/local/Cellar/ginac/1.7.5/include /usr/local/Cellar/cln/1.3.4/include /usr/local/Cellar/glpk/4.65/include)
include_directories(include include/STL parser/include)
file(GLOB_RECURSE SOURCES src/*.cpp src/STL/*.cpp parser/src/*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)

add_executable(sapo ${BISON_input_parser_OUTPUT_SOURCE} ${FLEX_input_scanner_OUTPUTS} ${SOURCES})

set(CMAKE_CXX_FLAGS "-O2")

target_compile_features(sapo PRIVATE cxx_range_for)
target_link_libraries(sapo ${PROJECT_LINK_LIBS} )
