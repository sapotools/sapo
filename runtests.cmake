macro(SAPO_TEST TEST_NAME TEST_DIR)
    execute_process(COMMAND bin/sapo ${TEST_DIR}/${TEST_NAME}.sil
                    RESULT_VARIABLE CMD_RESULT
                    OUTPUT_FILE ${TEST_NAME}.out)
    if(CMD_RESULT)
        message(FATAL_ERROR "Error on test ${TEST_NAME}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${TEST_DIR}/${TEST_NAME}.out ${TEST_NAME}.out
                    RESULT_VARIABLE CMD_RESULT)
    if(CMD_RESULT)
        message(FATAL_ERROR "Error on test ${TEST_NAME}")
    else()
        file(REMOVE ${TEST_NAME}.out)
    endif()
endmacro()

sapo_test(${TEST_NAME} ${TEST_DIR})