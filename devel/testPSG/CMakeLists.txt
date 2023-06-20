enable_testing()
add_test(
  NAME testPhyssurfShaoqi 
  COMMAND ${TARGET_NAME} data/main.dat
  WORKING_DIRECTORY ${INSTALL_PATH}
  )

