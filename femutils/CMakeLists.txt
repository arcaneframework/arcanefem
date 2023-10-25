cmake_minimum_required(VERSION 3.21)
project(FemTest1 LANGUAGES C CXX)

find_package(Arcane REQUIRED)

option(ENABLE_DEBUG_MATRIX "Enable Debug matrix instead of a sparse one" OFF)

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

target_compile_definitions(FemUtils PRIVATE $<$<BOOL:${ENABLE_DEBUG_MATRIX}>:ENABLE_DEBUG_MATRIX>)

target_include_directories(FemUtils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(FemUtils PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

arcane_add_arcane_libraries_to_target(FemUtils)

target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph)

set(FEMUTILS_HAS_PARALLEL_SOLVER FALSE)
set(FEMUTILS_HAS_PARALLEL_SOLVER_TRILINOS FALSE)
set(FEMUTILS_HAS_PARALLEL_SOLVER_HYPRE FALSE)
set(FEMUTILS_HAS_PARALLEL_SOLVER_PETSC FALSE)

if (TARGET Arcane::arcane_aleph_trilinos)
  set(FEMUTILS_HAS_PARALLEL_SOLVER TRUE)
  target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph_trilinos)
  message(STATUS "Trilinos backend is available")
  set(FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS TRUE)
endif()

if (TARGET Arcane::arcane_aleph_petsc)
  set(FEMUTILS_HAS_PARALLEL_SOLVER TRUE)
  target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph_petsc)
  message(STATUS "PETSc backend is available")
  set(FEMUTILS_HAS_SOLVER_BACKEND_PETSC TRUE)
endif()

if (TARGET Arcane::arcane_aleph_hypre)
  set(FEMUTILS_HAS_PARALLEL_SOLVER TRUE)
  target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph_hypre)
  message(STATUS "HYPRE backend is available")
  set(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE TRUE)
endif()

if(NOT FEMUTILS_HAS_PARALLEL_SOLVER)
  message(STATUS "WARNING: No Arcane parallel solver is available. Parallel linear solver system is disabled")
endif()

# Add a variable in cache to check which parallel solvers are availables
set(FEMUTILS_HAS_PARALLEL_SOLVER ${FEMUTILS_HAS_PARALLEL_SOLVER} CACHE BOOL "Is Parallel solver available" FORCE)
set(FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS ${FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS} CACHE BOOL "Is Trilinos solver available" FORCE)
set(FEMUTILS_HAS_SOLVER_BACKEND_PETSC ${FEMUTILS_HAS_SOLVER_BACKEND_PETSC} CACHE BOOL "Is PETSc solver available" FORCE)
set(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE ${FEMUTILS_HAS_SOLVER_BACKEND_HYPRE} CACHE BOOL "Is Hypre solver available" FORCE)
