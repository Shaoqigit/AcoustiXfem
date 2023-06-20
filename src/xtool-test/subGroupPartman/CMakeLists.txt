enable_testing()
add_test(
  NAME subGroupPartman 
  COMMAND  ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${MPIEXEC_PREFLAGS} ${INSTALL_PATH}/${TARGET_NAME} ${MPIEXEC_POSTFLAGS} data/main.dat
  WORKING_DIRECTORY ${INSTALL_PATH}
  )
add_test( 
  NAME ndiff_subGroupPartman 
  COMMAND  ${TESTNDIFF}  
  WORKING_DIRECTORY ${INSTALL_PATH} 
  )
set_tests_properties(ndiff_subGroupPartman PROPERTIES DEPENDS subGroupPartman  ) 
