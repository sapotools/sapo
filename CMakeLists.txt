cmake_minimum_required(VERSION 2.8.9)
project (sapo CXX)
 
#For shared libraries:
set ( PROJECT_LINK_LIBS libginac.so.6 libglpk.so)
link_directories( usr/local/lib )
 
include_directories(include include/models)
file(GLOB SOURCES ./src/*.cpp ./src/models/*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)
 
add_executable(sapo ${SOURCES})
target_compile_features(sapo PRIVATE cxx_range_for)
target_link_libraries(sapo ${PROJECT_LINK_LIBS} )


