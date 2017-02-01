cmake_minimum_required(VERSION 2.8.9)
project (sapo CXX)

#For shared libraries:
set ( PROJECT_LINK_LIBS ginac glpk)
#link_directories( usr/local/lib )

include_directories(include include/models include/STL)
file(GLOB_RECURSE SOURCES src/*.cpp src/models/*.cpp src/STL/*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ./bin)

add_executable(sapo ${SOURCES})

set(CMAKE_CXX_FLAGS "-O2")

target_compile_features(sapo PRIVATE cxx_range_for)
target_link_libraries(sapo ${PROJECT_LINK_LIBS} )
