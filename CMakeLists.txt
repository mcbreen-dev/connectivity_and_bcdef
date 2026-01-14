cmake_minimum_required(VERSION 3.20)

project(Connectivity_and_bcdef LANGUAGES C CXX)

if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Use the Homebrew CGNS install.
find_path(CGNS_INCLUDE_DIR cgnslib.h PATHS /opt/homebrew/include)
find_library(CGNS_LIBRARY cgns PATHS /opt/homebrew/lib)

if(NOT CGNS_INCLUDE_DIR OR NOT CGNS_LIBRARY)
  message(FATAL_ERROR "CGNS not found. Install via Homebrew or set CGNS_INCLUDE_DIR/CGNS_LIBRARY.")
endif()

add_library(CGNS::cgns-static UNKNOWN IMPORTED)
set_target_properties(CGNS::cgns-static PROPERTIES
  IMPORTED_LOCATION "${CGNS_LIBRARY}"
  INTERFACE_INCLUDE_DIRECTORIES "${CGNS_INCLUDE_DIR}"
)

add_subdirectory(BC_define)
