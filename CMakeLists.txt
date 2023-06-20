# ------------------------------------------------------------------------
# CMakeLists.txt for physsurfGeneric
# more documentation in  http://www.cmake.org/cmake/help/cmake2.6docs.html
# ------------------------------------------------------------------------
cmake_minimum_required(VERSION 2.6)
project(helmholtzNd)





message(STATUS "---------------------------")
message(STATUS "PROJECT: helmholtzNd --")
message(STATUS "---------------------------")

set(DEVROOT ${CMAKE_CURRENT_SOURCE_DIR}/../..)
list(APPEND CMAKE_MODULE_PATH  ${DEVROOT}/Xfiles/Util/cmakeUtil/  )
include(common_functions)
list(APPEND CMAKE_MODULE_PATH  ${DEVROOT} )
#include(ExternalLibrariesPaths)


#=======================================================================================
# libs and executables files root directory:
# variables CMAKE_INSTALL_DIR, CMAKE_SYSTEM_PREFIX_PATH and DESTDIR are not used.
#=======================================================================================
#set(LIBRARIES_INSTALL_PATH   "${DEVROOT}/LibDir"   CACHE INTERNAL "librairies root destination")
#set(EXECUTABLES_INSTALL_PATH "${DEVROOT}/ExecDir"  CACHE INTERNAL "executables root destination")

#BUILD_ROOT defini dans le CMakeLists principal
set(EXECUTABLES_INSTALL_PATH "${BUILD_ROOT}/ExecDir"  CACHE INTERNAL "executables root destination")


# Set C++ standard to cxx14 (because of bessel)
set(CMAKE_CXX_FLAGS "-std=c++14 -pthread" CACHE STRING "inforced c++ standard")

#=======================================================================================
# definition of options
#=======================================================================================
option(VERBOSE           "Verbose mode for CMake output"        ON )
# USE options are mandatory to find interface headers
option(USE_MUMPS         "use MUMPS interface"                  ON)

# options for solver : only one is on
option(MUMPSAPPLI        "APPLI using MUMPS  "                  OFF)
#
# to obtain compressed result files
option(GZIPOUTPUT "Use of boost iostream to compress output" OFF)

# use of hash map or hash set for map, set, ...
option(USE_HASH_TABLE "use hash table (default is ON)" ON)

# for special environement
option(MYENV "use special target for test and include in project attached project " OFF)


#option(HE_USE_OCTREE                     "Use Octree"                    ON )
option(HE_USE_EIGEN                      "Use Eigen"                     ON )


#=======================================================================================
# treatment of options
#=======================================================================================
# # setting cxx_compiler following chosen option
# set_cxx_compiler(EXTERNAL_INCLUDES EXTERNAL_LIBRARIES)
# # list of external package if needed
# set_external_def()
# set_internal_def()



#=======================================================================================
# Definitions
#=======================================================================================



if(NOT USE_XLEGACYSIMPLECUT)
   message(FATAL_ERROR AUTHOR_WARNING "
**************************************************************************
physsurfGeneric library need xLegacySimpleCut. You must set USE_XLEGACYSIMPLECUT to ON and (re)build xFEM with this option. If you do not need xPhysics, set BUILD_XPHYSICS to OFF
**************************************************************************
")
else()
  add_definitions( -DWITH_XLEGACYSIMPLECUT)
endif()

#=======================================================================================
# EXCEPTIONS
#=======================================================================================

#=======================================================================================
#=======================================================================================


# ---------------------------------------
# definition of the library 
# ---------------------------------------
file(GLOB src_files
  src/*.cc 
  src/*.h 
  )

# librairy can be SHARED or STATIC
set(BUILD_SHARED_LIBS "SHARED")
define_archos_suffixe(ARCHOS)

add_library(helmholtzNd  ${BUILD_SHARED_LIBS}  ${src_files}  )
set(LIBRARY_OUTPUT_PATH ${DEVROOT}/lib/${ARCHOS} ) 

find_package(xTool         REQUIRED)
find_package(xSolverBase            REQUIRED)
find_package(xInterfaceMumps      REQUIRED)
find_package(xInterfaceAOMDGeneral  REQUIRED)
find_package(Trellis      REQUIRED)
find_package(xFEM                   REQUIRED)
#find_package(xOctree                REQUIRED)
find_package(BLAS)
find_package(LAPACK2)
find_and_set(xInterfaceLapack xInterfaceLapack HE_EXTERNAL_INCLUDES HE_EXTERNAL_LIBRARIES)
find_package(METIS)
find_package(Boost                  REQUIRED)
find_package(xLegacySimpleCut       REQUIRED)
find_package(Sisl                  REQUIRED)
find_package(xGeomTools REQUIRED)
find_package(xQuadrature REQUIRED)
#find_package(xMeshTool REQUIRED)
#find_package(xInterfaceAOMDGeneral  REQUIRED)



find_and_set(xDistMesh xDistMesh HE_EXTERNAL_INCLUDES HE_EXTERNAL_LIBRARIES)
find_and_set(xInterfaceParMetis            xInterfaceParMetis            HE_EXTERNAL_INCLUDES  HE_EXTERNAL_LIBRARIES )

#find_and_set(xPhyssurfGeneric            xPhyssurfGeneric            HE_EXTERNAL_INCLUDES  HE_EXTERNAL_LIBRARIES )
find_package(xPhyssurfGeneric REQUIRED)
message("========================================PSG : ${xPSG_LIBRARIES}")
message("========================================PSG : ${xPSG_INCLUDE_DIR}")


if(HE_USE_OCTREE)
find_package(xOctree      REQUIRED)
list(APPEND HE_EXTERNAL_INCLUDES ${OCTREE_INCLUDE_DIR})
list(APPEND HE_EXTERNAL_LIBRARIES ${OCTREE_LIBRARIES})

find_package(xInterfacexOctree      REQUIRED)
list(APPEND HE_EXTERNAL_INCLUDES ${XoctreeInterface_INCLUDE_DIR})
list(APPEND HE_EXTERNAL_LIBRARIES ${XoctreeInterface_LIBRARIES})
endif()

if(HE_USE_EIGEN)
find_package(Eigen3)
list(APPEND HE_EXTERNAL_INCLUDES  ${EIGEN3_INCLUDE_DIR})
endif()



if (USE_MUMPS)
    find_package(xInterfaceMumps)
    list(APPEND  HE_EXTERNAL_INCLUDES     ${xInterfaceMumps_INCLUDE_DIR}    )
    add_definitions(-DUSE_MUMPS)# #define USE_MUMPS can be used in the source code
endif()

list(APPEND HE_EXTERNAL_INCLUDES  ${PARMETIS_INCLUDE_DIR} )
list(APPEND HE_EXTERNAL_INCLUDES     ${METIS_INCLUDE_DIR} )
list(APPEND HE_EXTERNAL_LIBRARIES ${PARMETIS_LIBRARIES}   )
list(APPEND HE_EXTERNAL_LIBRARIES    ${METIS_LIBRARIES}   )

list(APPEND HE_EXTERNAL_INCLUDES ${MPI_CXX_INCLUDE_PATH} )
list(APPEND HE_EXTERNAL_LIBRARIES ${MPI_CXX_LIBRARIES} )



list(APPEND HE_EXTERNAL_INCLUDES  
/develop/complex_bessel/install/include 
  ${xTool_INCLUDE_DIR} 
  ${xFEM_INCLUDE_DIR} 
  ${Trellis_INCLUDE_DIR} 
  ${xOctree_INCLUDE_DIR} 
  ${xInterfaceAOMDGeneral_INCLUDE_DIR} 
  ${xInterfacexOctree_INCLUDE_DIR} 
  ${xSolverBase_INCLUDE_DIR}  
  #${InterfaceMumps_INCLUDE_DIR}
  ${Lapack_INCLUDE_DIR}
  #${LIBIGL_INCLUDE_DIRS}
  ${xLegacySimpleCut_INCLUDE_DIR}
  ${xInterfaceLapack_INCLUDE_DIR}  
  ${LAPACK_INCLUDE_DIR} 
  ${xGeomTools_INCLUDE_DIR}
  ${xQuadrature_INCLUDE_DIR}
  ${LIBSISL_INCLUDE_DIR}
  ${xPSG_INCLUDE_DIR}

)  




include_directories( ${HE_EXTERNAL_INCLUDES}
)

set_target_properties(helmholtzNd PROPERTIES COMPILE_FLAGS "-w -Wno-deprecated ${CXXFLAGS}"  )


# ---------------------------------------
# adding devel/test
# ---------------------------------------
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/devel)

