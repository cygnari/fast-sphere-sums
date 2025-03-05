function (printvar var)
    message("${var}: ${${var}}")
endfunction()

cmake_minimum_required(VERSION 3.10)
project(FastSphereSums VERSION 1.0 LANGUAGES CXX C)
set(default_build_type Release)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}")
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)
option(BUILD_SHARED_LIBS "Build shared libraries, rather than static" ON)

if (BUILD_SHARED_LIBS)
  set(CMAKE_POSITION_INDEPENDENT_CODE ON)
endif()

set(CMAKE_CXX_FLAGS_RELEASE -O3)

# turn off deprecated env variable handling:
# cmake will now use $<Package>_ROOT env variables in its search for Package
# dependencies
CMAKE_POLICY(SET CMP0074 NEW)

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.24)
  # turn off deprecated cmake download timestamp handling
  # cmake will now set the local timestamp of a dependency to match its package
  # creation timestamp, rather than its download timestamp.
  CMAKE_POLICY(SET CMP0135 NEW)
endif()

find_package(MPI REQUIRED COMPONENTS C)
find_package(LAPACK REQUIRED)

add_subdirectory(src)
add_subdirectory(executables)

# CMake auto-generated config settings will go in this file
configure_file(fast-sphere-sums-config.h.in fast-sphere-sums-config.h)

# The following commands prepare FastSphereSums to be packaged as a library and used
# by other codes as a dependency.
#include(GNUInstallDirs)
#include(CMakePackageConfigHelpers)
#
# generate config file for other CMake projects to find this project
#configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
#  ${CMAKE_CURRENT_BINARY_DIR}/FastSphereSumsConfig.cmake
#  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/FastSphereSums"
#  NO_SET_AND_CHECK_MACRO
#  NO_CHECK_REQUIRED_COMPONENTS_MACRO
#)
# include version info so that client codes can ensure compatibility
#write_basic_package_version_file(
#  "${CMAKE_CURRENT_BINARY_DIR}/FastSphereSumsConfigVersion.cmake"
#  VERSION "${FastSphereSums_VERSION_MAJOR}.${FastSphereSums_VERSION_MINOR}"
#  COMPATIBILITY AnyNewerVersion
#)
# install configured package files
#install(FILES
#  ${CMAKE_CURRENT_BINARY_DIR}/FastSphereSumsConfig.cmake
#  ${CMAKE_CURRENT_BINARY_DIR}/FastSphereSumsConfigVersion.cmake
#  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FastSphereSums
#)
#install(EXPORT FastSphereSumsTargets
# FILE FastSphereSumsTargets.cmake
# DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/FastSphereSums
#)
# export configured package files to the build directory, too; this allows
# clients to use FastSphereSums from a FastSphereSums build directory (without having to
# install it).
#export(EXPORT FastSphereSumsTargets
# FILE "${CMAKE_CURRENT_BINARY_DIR}/FastSphereSumsTargets.cmake"
#)
