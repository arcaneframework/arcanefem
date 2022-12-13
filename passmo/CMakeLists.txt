cmake_minimum_required(VERSION 3.21)

project(Passmo LANGUAGES C CXX)
find_package(Arcane REQUIRED)

add_executable(Passmo
  AlephFemLinearSystemDoF.cc
  TypesElastodynamic.h
  Elastodynamic_axl.h
  ElastodynamicModule.h
  ElastodynamicModule.cc
  main.cc
  FemLinearSystemDoF.h
  FemLinearSystemDoF.cc
  FemUtils.h
  FemUtils.cc
  Integer3std.h
  utilFEM.h
  utilFEM.cc
  NodeDoFs_axl.h
  NodeDoFs.cc)

arcane_generate_axl(Elastodynamic)
arcane_generate_axl(NodeDoFs)
arcane_add_arcane_libraries_to_target(Passmo)

target_include_directories(Passmo PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Passmo.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Passmo.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Passmo PRIVATE Arcane::arcane_aleph)

if (TARGET Arcane::arcane_aleph_hypre)
  set(PASSMO_HAS_PARALLEL_SOLVER TRUE)
  target_link_libraries(Passmo PRIVATE Arcane::arcane_aleph_hypre)
else()
  message(STATUS "WARNING: Target Arcane::arcane_aleph_hypre is not available. Parallel linear solver system is disabled")
endif()


