find_package(Boost COMPONENTS unit_test_framework)

if(${Boost_unit_test_framework_FOUND})
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIRS})

    set(LIBSAPO_TESTS linear_algebra)
    foreach(TEST ${LIBSAPO_TESTS})
        ADD_EXECUTABLE( ${TEST} ${TEST}.cpp )
        target_link_libraries(${TEST} Sapo ${Boost_LIBRARIES})
        ADD_TEST(${TEST} ${TEST})
    endforeach()
endif()

