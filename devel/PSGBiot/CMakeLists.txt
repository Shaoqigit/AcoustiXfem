enable_testing()
add_test(
  NAME PSGBiot 
  COMMAND ${TARGET_NAME} 
  WORKING_DIRECTORY ${INSTALL_PATH}
  )

