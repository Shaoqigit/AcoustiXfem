enable_testing()
add_test(
  NAME runHelmholtz
  COMMAND ${TARGET_NAME} 
  WORKING_DIRECTORY ${INSTALL_PATH}
  )

