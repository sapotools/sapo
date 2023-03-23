cmake_minimum_required(VERSION 3.7)

project (sapo CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Werror ${CMAKE_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

set(THREADED_VERSION TRUE CACHE BOOL "Enable/disable threaded version compilation")

if(${THREADED_VERSION})
find_package(Threads)

if(${Threads_FOUND})
set(CMAKE_THREAD_PREFER_PTHREAD ON)

if(${CMAKE_VERSION} VERSION_LESS "3.12.0") 
add_definitions(-DWITH_THREADS)
else()
add_compile_definitions(WITH_THREADS)
endif()

set(PROJECT_LINK_LIBS ${PROJECT_LINK_LIBS} Threads::Threads)
endif()
else()
message("Threaded version disabled.")
endif()

execute_process(
        COMMAND git branch --show-current
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
   )

string(REGEX MATCH "^v[0-9]+[.][0-9]+.*$" SAPO_MINOR "${GIT_BRANCH}")

set(SAPO_VERSION "")

if(NOT SAPO_MINOR STREQUAL "")
   execute_process(
         COMMAND git describe --abbrev=0 --tags
         WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
         OUTPUT_VARIABLE GIT_TAG
         OUTPUT_STRIP_TRAILING_WHITESPACE
      )
   
   string(REGEX MATCH "^v[0-9]+[.][0-9]+[.][0-9]+$" SAPO_VERSION "${GIT_TAG}")
endif()

if(SAPO_VERSION STREQUAL "")
   execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}
        OUTPUT_VARIABLE SAPO_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
   )

   set(SAPO_VERSION "commit ${SAPO_VERSION}")
endif()

message("This is Sapo ${SAPO_VERSION}")

set(RUNTEST_FILE "${CMAKE_SOURCE_DIR}/scripts/run-tests.bash" )

file(READ "${RUNTEST_FILE}.orig" FILEDATA)
string(REGEX REPLACE "CMAKE_MAKE_PROGRAM" "${CMAKE_MAKE_PROGRAM}" FILEDATA "${FILEDATA}")
file(WRITE ${RUNTEST_FILE} "${FILEDATA}")

execute_process(COMMAND sh ${CMAKE_SOURCE_DIR}/scripts/install-hooks.bash )

find_program( CLANG_FORMAT_EXE NAMES "clang-format")
if(NOT CLANG_FORMAT_EXE)
   message(STATUS "clang-format not found.")
else()
   message(STATUS "clang-format found: ${CLANG_FORMAT_EXE}")
endif()

enable_testing()

add_subdirectory(libSapo)
add_subdirectory(sapo)
