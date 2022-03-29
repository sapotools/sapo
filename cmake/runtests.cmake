macro(SAPO_TEST TEST_NAME TEST_DIR)
    execute_process(COMMAND bin/sapo -t ${TEST_DIR}/${TEST_NAME}.sil
                    RESULT_VARIABLE CMD_RESULT
                    OUTPUT_FILE ${TEST_NAME}.out
	    	    ERROR_FILE ${TEST_NAME}.err)

    if(CMD_RESULT)
    	execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${TEST_DIR}/${TEST_NAME}.err ${TEST_NAME}.err
                        RESULT_VARIABLE CMD_RESULT)

    	if(CMD_RESULT)
            message(FATAL_ERROR "Error on test ${TEST_NAME}")
        else()
            file(REMOVE ${TEST_NAME}.err)
            file(REMOVE ${TEST_NAME}.out)
        endif()
    else()
       	file(REMOVE ${TEST_NAME}.err)

    	execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files ${TEST_DIR}/${TEST_NAME}.out ${TEST_NAME}.out
                        RESULT_VARIABLE CMD_RESULT)
    	if(CMD_RESULT)
            message(FATAL_ERROR "Error on test ${TEST_NAME}")
    	else()
            file(REMOVE ${TEST_NAME}.out)
    	endif()
    endif()
endmacro()

sapo_test(${TEST_NAME} ${TEST_DIR})
