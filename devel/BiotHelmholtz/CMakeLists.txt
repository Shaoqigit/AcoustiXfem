enable_testing()
add_test(
  NAME BiotHelmholtz
  COMMAND ${TARGET_NAME} 
  WORKING_DIRECTORY ${INSTALL_PATH}
  )

