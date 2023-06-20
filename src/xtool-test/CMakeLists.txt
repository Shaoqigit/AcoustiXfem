cmake_minimum_required(VERSION 2.6)
set(text       "-------------------------- --")
message_verbose(text)  
message(STATUS "building Xtest/xtools_test --")
set(text       "-------------------------- --")
message_verbose(text)  
 

set_external_def()

 
find_and_set(xTool                         xTool                         EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(Boost                         Boost                         EXTERNAL_INCLUDES  NOTHING             REQUIRED )
find_and_set(xDistMesh                     xDistMesh                     EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xTemplateRefineMesh           xTemplateRefineMesh           EXTERNAL_INCLUDES  NOTHING             REQUIRED )
find_and_set(xInterfacexTemplateRefineMesh xInterfacexTemplateRefineMesh EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xInterfaceAOMDGeneral         xInterfaceAOMDGeneral         EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(XtestCommon                   XTESTCOMMON                   EXTERNAL_INCLUDES  NOTHING             REQUIRED )
find_and_set(Trellis                       Trellis                       EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xExport                       xExport                       EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xTensor                       xTensor                       EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xFEM                          xFEM                          EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xQuadrature                   xQuadrature                   EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xMapping                      xMapping                      EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xSolverBase                   xSolverBase                   EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xLegacySimpleCut              xLegacySimpleCut              EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xInterfaceLapack              xInterfaceLapack              EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(xInterfaceParMetis            xInterfaceParMetis            EXTERNAL_INCLUDES  EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(LAPACK                        LAPACK                        NOTHING            EXTERNAL_LIBRARIES  REQUIRED )
find_and_set(BLAS2                         BLAS                          NOTHING            EXTERNAL_LIBRARIES  REQUIRED )
list(APPEND EXTERNAL_INCLUDES  ${MPI_CXX_INCLUDE_PATH} )
list(APPEND EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES} )

find_package(METIS       REQUIRED )
set_from_prefix(PARMETIS EXTERNAL_INCLUDES EXTERNAL_LIBRARIES )
set_from_prefix(METIS    EXTERNAL_INCLUDES EXTERNAL_LIBRARIES )

enable_testing()


set(LIST
    ${CMAKE_CURRENT_SOURCE_DIR}/subGroupPartman
     )


create_tests_from_list(${LIST})
