cmake_minimum_required(VERSION 3.20)

project(Connectivity_and_bcdef LANGUAGES C CXX)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

option(AUTO_DEPS "Download/build CGNS + HDF5 if not found" ON)
set(CGNS_VERSION "4.5.0" CACHE STRING "CGNS version for AUTO_DEPS")
set(HDF5_VERSION "1.14.3" CACHE STRING "HDF5 version for AUTO_DEPS")

if(APPLE)
  include(${CMAKE_CURRENT_LIST_DIR}/macos/cmake/Defaults.cmake)
elseif(UNIX)
  include(${CMAKE_CURRENT_LIST_DIR}/linux/cmake/Defaults.cmake)
endif()

set(CGNS_TARGET "")

# Try to load HDF5 config if it exists (needed by CGNS config targets).
find_package(HDF5 CONFIG QUIET)

# Prefer a CMake package if available.
find_package(cgns CONFIG QUIET)
if(TARGET CGNS::cgns_shared)
  set(CGNS_TARGET CGNS::cgns_shared)
elseif(TARGET CGNS::cgns_static)
  set(CGNS_TARGET CGNS::cgns_static)
endif()

# Fallback to explicit CGNS_INCLUDE_DIR/CGNS_LIBRARY if provided.
if(NOT CGNS_TARGET AND CGNS_INCLUDE_DIR AND CGNS_LIBRARY)
  add_library(CGNS::cgns_static UNKNOWN IMPORTED)
  set_target_properties(CGNS::cgns_static PROPERTIES
    IMPORTED_LOCATION "${CGNS_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${CGNS_INCLUDE_DIR}"
  )
  set(CGNS_TARGET CGNS::cgns_static)
endif()

# Auto-download dependencies if requested.
if(NOT CGNS_TARGET AND AUTO_DEPS)
  include(FetchContent)
  message(STATUS "CGNS not found; fetching CGNS and HDF5")

  set(HDF5_BUILD_CPP_LIB OFF CACHE BOOL "" FORCE)
  set(HDF5_BUILD_FORTRAN OFF CACHE BOOL "" FORCE)
  set(HDF5_BUILD_TOOLS OFF CACHE BOOL "" FORCE)
  set(HDF5_BUILD_UTILS OFF CACHE BOOL "" FORCE)
  set(HDF5_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(HDF5_BUILD_HL_LIB OFF CACHE BOOL "" FORCE)
  set(HDF5_ENABLE_Z_LIB_SUPPORT OFF CACHE BOOL "" FORCE)
  set(HDF5_ENABLE_SZIP_SUPPORT OFF CACHE BOOL "" FORCE)

  FetchContent_Declare(hdf5
    URL https://github.com/HDFGroup/hdf5/archive/refs/tags/hdf5-${HDF5_VERSION}.tar.gz
  )
  FetchContent_MakeAvailable(hdf5)

  set(HDF5_ROOT "${hdf5_BINARY_DIR}" CACHE PATH "" FORCE)
  set(HDF5_DIR "${hdf5_BINARY_DIR}" CACHE PATH "" FORCE)

  set(CGNS_ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
  set(CGNS_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
  set(CGNS_BUILD_SHARED ON CACHE BOOL "" FORCE)
  set(CGNS_USE_HDF5 ON CACHE BOOL "" FORCE)
  set(CGNS_ENABLE_HDF5 ON CACHE BOOL "" FORCE)

  FetchContent_Declare(cgns
    URL https://github.com/CGNS/CGNS/archive/refs/tags/v${CGNS_VERSION}.tar.gz
  )
  FetchContent_MakeAvailable(cgns)

  if(TARGET CGNS::cgns_shared)
    set(CGNS_TARGET CGNS::cgns_shared)
  elseif(TARGET CGNS::cgns_static)
    set(CGNS_TARGET CGNS::cgns_static)
  elseif(TARGET cgns_shared)
    add_library(CGNS::cgns_shared ALIAS cgns_shared)
    set(CGNS_TARGET CGNS::cgns_shared)
  elseif(TARGET cgns_static)
    add_library(CGNS::cgns_static ALIAS cgns_static)
    set(CGNS_TARGET CGNS::cgns_static)
  elseif(TARGET cgns)
    add_library(CGNS::cgns_shared ALIAS cgns)
    set(CGNS_TARGET CGNS::cgns_shared)
  endif()
endif()

if(NOT CGNS_TARGET)
  message(FATAL_ERROR "CGNS not found. Install CGNS (with HDF5) or configure with AUTO_DEPS=ON.")
endif()

set(CGNS_TARGET "${CGNS_TARGET}" CACHE STRING "CGNS CMake target to link" FORCE)

add_subdirectory(BC_define)
