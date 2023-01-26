find_package(Boost COMPONENTS unit_test_framework)

if(${Boost_unit_test_framework_FOUND})
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})

    set(LIBSAPO_TESTS linear_algebra symbolic_algebra linear_systems
                    polytopes parallelotopes bundles dynamical_systems
                    evolver sets_unions sticky_unions Bernstein
                    simplex)
    foreach(TEST ${LIBSAPO_TESTS})
        ADD_EXECUTABLE( test_${TEST} ${TEST}.cpp )
        target_link_libraries(test_${TEST} Sapo ${Boost_LIBRARIES})
        ADD_TEST(${TEST} test_${TEST})
    endforeach()
endif()

