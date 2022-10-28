CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT( social_octomap )

include(CTest)
include(GNUInstallDirs)

# version (e.g. for packaging)
set(SOCIAL_OCTOMAP_MAJOR_VERSION 1)
set(SOCIAL_OCTOMAP_MINOR_VERSION 9)
set(SOCIAL_OCTOMAP_PATCH_VERSION 8)
set(SOCIAL_OCTOMAP_VERSION ${SOCIAL_OCTOMAP_MAJOR_VERSION}.${SOCIAL_OCTOMAP_MINOR_VERSION}.${SOCIAL_OCTOMAP_PATCH_VERSION})
set(SOCIAL_OCTOMAP_SOVERSION ${SOCIAL_OCTOMAP_MAJOR_VERSION}.${SOCIAL_OCTOMAP_MINOR_VERSION})
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
  if(POLICY CMP0042)
    # Enable MACOSX_RPATH by default.
    cmake_policy(SET CMP0042 NEW)
  endif(POLICY CMP0042)
endif(COMMAND cmake_policy)

SET (CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/CMakeModules")

# COMPILER SETTINGS (default: Release) and flags
INCLUDE(CompilerSettings)

# SOCIAL_OCTOMAP_OMP = enable OpenMP parallelization (experimental, defaults to OFF)
SET(SOCIAL_OCTOMAP_OMP FALSE CACHE BOOL "Enable/disable OpenMP parallelization")
IF(DEFINED ENV{SOCIAL_OCTOMAP_OMP})
  SET(SOCIAL_OCTOMAP_OMP $ENV{SOCIAL_OCTOMAP_OMP})
ENDIF(DEFINED ENV{SOCIAL_OCTOMAP_OMP})
IF(SOCIAL_OCTOMAP_OMP)
  FIND_PACKAGE( OpenMP REQUIRED)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ENDIF(SOCIAL_OCTOMAP_OMP)

# Set output directories for libraries and executables
SET( BASE_DIR ${CMAKE_SOURCE_DIR} )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BASE_DIR}/lib )
SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BASE_DIR}/lib )
SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BASE_DIR}/bin )
# output dirs for multi-config builds (MSVC)
foreach( OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES} )
  STRING( TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG )
  SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${BASE_DIR}/lib )
  SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${BASE_DIR}/lib )
  SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG} ${BASE_DIR}/bin )
endforeach( OUTPUTCONFIG CMAKE_CONFIGURATION_TYPES )

set(INCLUDE_DIRS "${PROJECT_SOURCE_DIR}/include")
INCLUDE_DIRECTORIES(${INCLUDE_DIRS})

LINK_DIRECTORIES(${CMAKE_LIBRARY_OUTPUT_DIRECTORY})

# Installation

set(INSTALL_TARGETS_DEFAULT_ARGS
	RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
	LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
	ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)

ADD_SUBDIRECTORY( src/math )
ADD_SUBDIRECTORY( src )

file(GLOB social_octomap_HDRS ${PROJECT_SOURCE_DIR}/include/social_octomap/*.h ${PROJECT_SOURCE_DIR}/include/social_octomap/*.hxx)
install(FILES ${social_octomap_HDRS} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/social_octomap")
file(GLOB social_octomap_math_HDRS ${PROJECT_SOURCE_DIR}/include/social_octomap/math/*.h)
install(FILES ${social_octomap_math_HDRS} DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/social_octomap/math")

# Install package.xml (catkin/ament/rosdep)
install(FILES package.xml DESTINATION "${CMAKE_INSTALL_DATADIR}/social_octomap")

# Allows Colcon to find non-Ament packages when using workspace underlays
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/share/ament_index/resource_index/packages/${PROJECT_NAME} "")
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/share/ament_index/resource_index/packages/${PROJECT_NAME} DESTINATION share/ament_index/resource_index/packages)

# uninstall target
configure_file(
    "${PROJECT_SOURCE_DIR}/CMakeModules/CMakeUninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)

# Export the package for use from the build-tree
# (this registers the build-tree with a global CMake-registry)
export(PACKAGE social_octomap)
 
# Create a social_octomap-config.cmake file for the use from the build tree
set(SOCIAL_OCTOMAP_INCLUDE_DIRS "${INCLUDE_DIRS}")
set(SOCIAL_OCTOMAP_LIB_DIR "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
# Set library names as absolute paths
# Windows, spec. MSVC requires the .lib suffix for imported libs
IF(WIN32)
  set(SOCIAL_OCTOMAP_LIBRARY
    "${CMAKE_IMPORT_LIBRARY_PREFIX}social_octomap${CMAKE_IMPORT_LIBRARY_SUFFIX}"
  )
  set(OCTOMATH_LIBRARY
    "${CMAKE_IMPORT_LIBRARY_PREFIX}octomath${CMAKE_IMPORT_LIBRARY_SUFFIX}"
  )
ELSE()
  set(SOCIAL_OCTOMAP_LIBRARY
    "${CMAKE_SHARED_LIBRARY_PREFIX}social_octomap${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )
  set(OCTOMATH_LIBRARY
    "${CMAKE_SHARED_LIBRARY_PREFIX}octomath${CMAKE_SHARED_LIBRARY_SUFFIX}"
  )
ENDIF()

# not used right now (export depends?)
#set(SOCIAL_OCTOMAP_CMAKE_DIR "${PROJECT_BINARY_DIR}")
include(CMakePackageConfigHelpers)

CONFIGURE_PACKAGE_CONFIG_FILE(
  social_octomap-config.cmake.in
  "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap/social_octomap-config.cmake"
  PATH_VARS SOCIAL_OCTOMAP_INCLUDE_DIRS SOCIAL_OCTOMAP_LIB_DIR
  INSTALL_DESTINATION "${CMAKE_INSTALL_FULL_DATADIR}/social_octomap")

WRITE_BASIC_PACKAGE_VERSION_FILE(
  "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/cmake/social_octomap/social_octomap-config-version.cmake"
  VERSION ${SOCIAL_OCTOMAP_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Create a social_octomap-config.cmake file for the use from the install tree
# and install it
set(SOCIAL_OCTOMAP_INCLUDE_DIRS "${CMAKE_INSTALL_FULL_INCLUDEDIR}")
set(SOCIAL_OCTOMAP_LIB_DIR "${CMAKE_INSTALL_FULL_LIBDIR}")
#set(SOCIAL_OCTOMAP_CMAKE_DIR "${INSTALL_DATA_DIR}/FooBar/CMake")

set(SOCIAL_OCTOMAP_INCLUDE_TARGETS
  "include(\${CMAKE_CURRENT_LIST_DIR}/social_octomap-targets.cmake)")

CONFIGURE_PACKAGE_CONFIG_FILE(
  social_octomap-config.cmake.in
  "${PROJECT_BINARY_DIR}/InstallFiles/social_octomap-config.cmake"
  PATH_VARS SOCIAL_OCTOMAP_INCLUDE_DIRS SOCIAL_OCTOMAP_LIB_DIR
  INSTALL_DESTINATION "${CMAKE_INSTALL_FULL_DATADIR}/social_octomap")

WRITE_BASIC_PACKAGE_VERSION_FILE(
  ${PROJECT_BINARY_DIR}/InstallFiles/social_octomap-config-version.cmake
  VERSION ${SOCIAL_OCTOMAP_VERSION}
  COMPATIBILITY AnyNewerVersion)

install(FILES
  "${PROJECT_BINARY_DIR}/InstallFiles/social_octomap-config.cmake"
  "${PROJECT_BINARY_DIR}/InstallFiles/social_octomap-config-version.cmake" 
  DESTINATION "${CMAKE_INSTALL_DATADIR}/social_octomap")

# Write pkgconfig-file:
include(InstallPkgConfigFile)
install_pkg_config_file(social_octomap
    CFLAGS
    LIBS -lsocial_octomap -loctomath
    REQUIRES
    VERSION ${SOCIAL_OCTOMAP_VERSION})

# Documentation 
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/social_octomap.dox.in ${CMAKE_CURRENT_BINARY_DIR}/social_octomap.dox @ONLY)
  ADD_CUSTOM_TARGET(docs ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/social_octomap.dox
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    COMMENT "Generating documentation (Doxygen) at location ${PROJECT_SOURCE_DIR}/doc/html/")
ENDIF(DOXYGEN_FOUND)

# Needs to be last statement:
INCLUDE(CPackSettings)

# Finished:
MESSAGE ("\n")
MESSAGE (STATUS "Compile social_octomap using: make")
MESSAGE (STATUS "Install social_octomap using: make install")
MESSAGE (STATUS "    (be sure to set the correct CMAKE_INSTALL_PREFIX before)")
MESSAGE (STATUS "Compile API-documentation using: make docs\n")
