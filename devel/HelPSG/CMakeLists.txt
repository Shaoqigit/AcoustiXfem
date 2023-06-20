enable_testing()
add_test(
  NAME HelPSG
  COMMAND ${TARGET_NAME} 
  WORKING_DIRECTORY ${INSTALL_PATH}
  )

