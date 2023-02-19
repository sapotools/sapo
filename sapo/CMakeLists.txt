include_directories(include include/)

find_package(FLEX 2.6.3)
find_package(BISON 3.6.0)

if (FLEX_FOUND)
if (BISON_FOUND)

set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" IS_SYSTEM_DIR)
if("${IS_SYSTEM_DIR}" STREQUAL "-1")
    set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif("${IS_SYSTEM_DIR}" STREQUAL "-1")

flex_target(input_scanner scanner.l
            ${CMAKE_CURRENT_SOURCE_DIR}/src/scanner.cpp)
bison_target(input_parser parser.yy
             ${CMAKE_CURRENT_SOURCE_DIR}/src/parser.cpp
             DEFINES_FILE ${CMAKE_CURRENT_SOURCE_DIR}/include/parser.h)

include_directories(${CMAKE_SOURCE_DIR}/libSapo/include ${CMAKE_SOURCE_DIR}/libSapo/include/STL)

file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

add_executable(sapo ${BISON_input_parser_OUTPUT_SOURCE}
	            ${FLEX_input_scanner_OUTPUTS} ${SOURCES})

target_link_libraries(sapo Sapo)

install(TARGETS sapo RUNTIME DESTINATION bin)

if(CLANG_FORMAT_EXE)
   file(GLOB HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/include/*.h)

   add_custom_target( format-sapo-code COMMAND ${CLANG_FORMAT_EXE} -i
                             -style=file ${SOURCES} ${HEADERS} )
endif()

set(sapo_tests Ebola Influenza Influenza-no-splits Influenza-splits
               LotkaVolterra Phosphorelay Quadcopter Rossler SIRp SIR
	            VanDerPol SIR_vax_assume incomplete-directions 
               invariant_validation invariant_validation_fail
               Mandelbrot_invariant rational_polynomial
               function_test)

add_test(sapo_compilation "${CMAKE_COMMAND}" --build ${CMAKE_BINARY_DIR}
	                   --target sapo -j 4 -- )
set(SAPO_EXEC "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/sapo")
foreach(test_name ${sapo_tests})
   add_test(NAME ${test_name}
	    COMMAND ${CMAKE_COMMAND} -DSAPO_EXEC=${SAPO_EXEC} 
          -DTEST_NAME=${test_name} -DTEST_DIR=${CMAKE_CURRENT_SOURCE_DIR}/examples
          -DTEST_RESULTS_DIR=${CMAKE_CURRENT_SOURCE_DIR}/tests
			 -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/runtests.cmake)
   set_tests_properties(${test_name} PROPERTIES FIXTURES_SETUP sapo_compilation)
endforeach()

else(BISON_FOUND)
message("Bison is not available: sapo standalone application will not be compiled")

endif(BISON_FOUND)
else(FLEX_FOUND)

message("Flex is not available: sapo standalone application will not be compiled")

endif(FLEX_FOUND)
