cmake_minimum_required(VERSION 3.21)
project(FemTest1 LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_library(FemUtils
  FemUtils.h
  FemUtils.cc
  NodeLinearSystem.h
  NodeLinearSystem.cc
  DoFLinearSystem.h
  DoFLinearSystem.cc
  FemDoFsOnNodes.h
  FemDoFsOnNodes.cc
  AlephNodeLinearSystem.cc
  AlephDoFLinearSystem.cc
  IDoFLinearSystemFactory.h
  AlephDoFLinearSystemFactory_axl.h
  SequentialBasicDoFLinearSystemFactory_axl.h
)
arcane_generate_axl(AlephDoFLinearSystemFactory)
arcane_generate_axl(SequentialBasicDoFLinearSystemFactory)

target_include_directories(FemUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FemUtils PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

arcane_add_arcane_libraries_to_target(FemUtils)

target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph)
if (TARGET Arcane::arcane_aleph_trilinos)
  target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph_trilinos)
endif()
if (TARGET Arcane::arcane_aleph_petsc)
  target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph_petsc)
endif()
if (TARGET Arcane::arcane_aleph_hypre)
  set(FEMUTILS_HAS_PARALLEL_SOLVER TRUE)
  target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph_hypre)
else()
  set(FEMUTILS_HAS_PARALLEL_SOLVER FALSE)
  message(STATUS "WARNING: Target Arcane::arcane_aleph_hypre is not available. Parallel linear solver system is disabled")
endif()

# Add a variable in cache to check if parallel solver is available
set(FEMUTILS_HAS_PARALLEL_SOLVER ${FEMUTILS_HAS_PARALLEL_SOLVER} CACHE BOOL "Is Parallel solver available" FORCE)
