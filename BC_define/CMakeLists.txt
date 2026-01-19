cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED CGNS_TARGET)
  message(FATAL_ERROR "CGNS_TARGET not set. Configure from the top-level CMakeLists.txt")
endif()

set(BC_DEFINE_COMMON_FLAGS -O3 -march=native -ffast-math)

function(bc_define_target_setup target_name)
  target_compile_features(${target_name} PRIVATE cxx_std_20)
  target_compile_options(${target_name} PRIVATE ${BC_DEFINE_COMMON_FLAGS})
  target_compile_definitions(${target_name} PRIVATE VERIFY_TRANSFORM)
  target_include_directories(${target_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
  target_link_libraries(${target_name} PRIVATE ${CGNS_TARGET})
  set_target_properties(${target_name} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin"
  )
endfunction()

set(BC_DEFINE_SOURCES
  src/bc_define_main.cpp
  src/bc_define_bcspec.cpp
  src/bc_define_cgns.cpp
  src/bc_define_plot3d.cpp
  src/bc_define_helpers.cpp
  src/mesh_io.cpp
  src/cgns_cleanup.cpp
  src/plot3d_io.cpp
  src/connectivity_core.cpp
  src/connectivity_face_build.cpp
  src/connectivity_hash.cpp
  src/connectivity_verify.cpp
  src/connectivity_utils.cpp
  src/connectivity_writer.cpp
  src/logger.cpp
)

add_executable(bc_define ${BC_DEFINE_SOURCES})
bc_define_target_setup(bc_define)

set(BC_DUMP_SOURCES
  src/bc_dump.cpp
  src/mesh_io.cpp
  src/cgns_cleanup.cpp
  src/logger.cpp
)

add_executable(bc_dump ${BC_DUMP_SOURCES})
bc_define_target_setup(bc_dump)
