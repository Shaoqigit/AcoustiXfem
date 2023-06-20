cmake_minimum_required(VERSION 2.6)
set(text "----------------------------")
message_verbose(text)  
message(STATUS "building helmholtzNd/devel")
set(text "----------------------------")
message_verbose(text)  
 
# ---------------------------------------
# setting options
# option can be ON or OFF
# the test if() is true if it's ON
# ---------------------------------------
# parallel
# # tau automatic profiling
# option(USE_TAU "Use Tau automatic profiling/trace compiler" OFF)
# # solver choice
# option(MUMPSAPPLI "APPLI using MUMPS  " ON)
# # use of gzip for output
# option(GZIPOUTPUT "Use of boost iostream to compress output" OFF)


# --------------------------------------------
# setting cxx_compiler following chosen option
# --------------------------------------------
# set_cxx_compiler(EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
#project(xfem_par_test)

# -------------------------
# treatment of some options
# -------------------------
# set_internal_def()
set_external_def()

# ---------------------------------------
# Find external package
# ---------------------------------------
find_package(xSolverBase            REQUIRED)
#find_package(xInterfaceAOMDGeneral            REQUIRED)
set_from_prefix(SolverBase EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
set_from_prefix(Lapack EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)

find_and_set(xTool xTool EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
find_and_set(xDistMesh DISTMESH EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
find_and_set(xInterfaceAOMDGeneral xInterfaceAOMDGeneral  EXTERNAL_INCLUDES EXTERNAL_LIBRARIES REQUIRED)
find_and_set(Trellis Trellis EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
find_and_set(xInterfaceLapack xInterfaceLapack EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
find_and_set(xFEM xFEM EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
find_and_set(xExport xExport EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
find_and_set(xDistMesh xDistMesh EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
set_from_prefix(xSolverBase EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)





if(HE_USE_OCTREE)
find_package(xOctree)
list(APPEND EXTERNAL_INCLUDES ${xOctree_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBRARIES ${xOctree_LIBRARIES})

find_package(xInterfacexOctree)
list(APPEND EXTERNAL_INCLUDES ${xInterfacexOctree_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBRARIES ${xInterfacexOctree_LIBRARIES})
endif()



#To compile (and especially link) MPI...
list(APPEND EXTERNAL_INCLUDES ${MPI_CXX_INCLUDE_PATH} )
list(APPEND EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES} )
#list(APPEND EXTERNAL_LIBRARIES ${CGAL_LIBRARIES} )


# For xReader...
list(APPEND EXTERNAL_INCLUDES ${xInterfaceAOMDGeneral_INCLUDE_DIR} )
list(APPEND EXTERNAL_LIBRARIES ${xInterfaceAOMDGeneral_LIBRARIES} )

if (USE_MUMPS)
  add_definitions( -DWITH_DMUMPS)
  find_package(xInterfaceMumps)
  set_from_prefix(xInterfaceMumps EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
  find_and_set(MUMPS MUMPS NOTHING EXTERNAL_LIBRARIES)
  list(APPEND EXTERNAL_LIBRARIES ${MPI_Fortran_LIBRARIES} )
endif()


find_and_set(LAPACK LAPACK NOTHING EXTERNAL_LIBRARIES)
find_and_set(BLAS BLAS NOTHING EXTERNAL_LIBRARIES)

find_package(METIS)
find_package(xLegacySimpleCut REQUIRED)
find_and_set(xMapping              xMapping      EXTERNAL_INCLUDES EXTERNAL_LIBRARIES  REQUIRED )
find_package(xPhyssurfGeneric REQUIRED)


list(APPEND EXTERNAL_INCLUDES    ${CMAKE_CURRENT_SOURCE_DIR}/../src 
${xInterfaceAOMDGeneral_INCLUDE_DIR})
list(APPEND EXTERNAL_LIBRARIES   helmholtzNd   )

list(APPEND EXTERNAL_INCLUDES     ${PARMETIS_INCLUDE_DIR} )
list(APPEND EXTERNAL_INCLUDES     ${METIS_INCLUDE_DIR} )
list(APPEND EXTERNAL_LIBRARIES    ${PARMETIS_LIBRARIES} )
list(APPEND EXTERNAL_LIBRARIES    ${METIS_LIBRARIES} )
list(APPEND EXTERNAL_LIBRARIES    ${xInterfaceLapack_LIBRARIES} )#Bizarre, ne fct pas avec find_and_set
list(APPEND EXTERNAL_LIBRARIES ${xLegacySimpleCut_LIBRARIES}  )
#list(APPEND EXTERNAL_LIBRARIES    -lboost_program_options )
#list(APPEND EXTERNAL_LIBRARIES    -lboost_system )
#list(APPEND EXTERNAL_LIBRARIES    -lboost_filesystem )
list(APPEND EXTERNAL_LIBRARIES    -lz )

list(APPEND EXTERNAL_LIBRARIES
    ${xTool_LIBRARIES} 
     ${xFEM_LIBRARIES} 
     ${xExport_LIBRARIES} 
  ${Trellis_LIBRARIES} 
  ${xInterfaceAOMDGeneral_LIBRARIES} 
  ${xSolverBase_LIBRARIES}  
  ${InterfaceMumps_LIBRARIES}
  ${InterfaceSuperLu_LIBRARIES}
  ${Lapack_LIBRARIES}
  #${xLegacySimpleCut_LIBRARIES}
    #${xInterfaceLapack_LIBRARIES}  
  ${LAPACK_LIBRARIES} 
    ${LIBSISL_LIBRARY}
    ${xQuadrature_LIBRARIES}
    ${BLAS_LIBRARIES}
    ${xPSG_LIBRARIES}
)
list(APPEND EXTERNAL_LIBRARIES  /develop/complex_bessel/install/lib/libcomplex_bessel.so )
enable_testing()


set(LIST
	#    ${CMAKE_CURRENT_SOURCE_DIR}/inutile
    #${CMAKE_CURRENT_SOURCE_DIR}/runHelmholtz
	#${CMAKE_CURRENT_SOURCE_DIR}/BiotHelmholtz
	${CMAKE_CURRENT_SOURCE_DIR}/HelPSG
	#${CMAKE_CURRENT_SOURCE_DIR}/PSGBiot
)

create_tests_from_list(${LIST})
