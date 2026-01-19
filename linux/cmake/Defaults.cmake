# Linux defaults for CGNS lookup
if(NOT CGNS_INCLUDE_DIR)
  find_path(CGNS_INCLUDE_DIR cgnslib.h PATHS /usr/local/include /usr/include)
endif()
if(NOT CGNS_LIBRARY)
  find_library(CGNS_LIBRARY cgns PATHS /usr/local/lib /usr/lib /usr/lib64)
endif()
