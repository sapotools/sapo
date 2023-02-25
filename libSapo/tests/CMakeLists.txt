find_package(Boost COMPONENTS unit_test_framework)

find_package(GMP)

if(${Boost_unit_test_framework_FOUND})
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})
    set(TEST_LIBRARIES Sapo ${Boost_LIBRARIES})

    if (${GMP_FOUND})
        set(TEST_LIBRARIES ${TEST_LIBRARIES} GMP::gmpxx GMP::gmp)
        include_directories(${GMPXX_INCLUDE_DIR} ${GMP_INCLUDE_DIR})

        if(${CMAKE_VERSION} VERSION_LESS "3.12.0") 
            add_definitions(-DHAVE_GMP)
        else()
            add_compile_definitions(HAVE_GMP)
        endif()
    endif()

    set(LIBSAPO_TESTS floating_points approximations linear_algebra
                    simplex linear_systems symbolic_algebra 
                    Bernstein polytopes parallelotopes bundles
                    evolver ode sets_unions sticky_unions
                    discrete_systems)
    foreach(TEST ${LIBSAPO_TESTS})
        ADD_EXECUTABLE( test_${TEST} ${TEST}.cpp )
        if (${GMP_FOUND})
            target_link_libraries(test_${TEST} ${TEST_LIBRARIES})
        else()
            target_link_libraries(test_${TEST} ${TEST_LIBRARIES})
        endif()
        ADD_TEST(${TEST} test_${TEST})
    endforeach()
endif()

