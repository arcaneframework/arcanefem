option(ENABLE_DEBUG_MATRIX "Enable Debug matrix instead of a sparse one" OFF)

add_library(FemUtils
  FemUtils.h
  FemUtils.cc
  DoFLinearSystem.h
  DoFLinearSystem.cc
  CooFormatMatrix.h
  CsrFormatMatrix.h
  CsrFormatMatrix.cc
  FemDoFsOnNodes.h
  FemDoFsOnNodes.cc
  AlephDoFLinearSystem.cc
  IDoFLinearSystemFactory.h
  AlephDoFLinearSystemFactory_axl.h
  SequentialBasicDoFLinearSystemFactory_axl.h
  HypreDoFLinearSystemFactory_axl.h
)

arcane_generate_axl(AlephDoFLinearSystemFactory)
arcane_generate_axl(SequentialBasicDoFLinearSystemFactory)
arcane_generate_axl(HypreDoFLinearSystemFactory)

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

if(0)
if (TARGET Arcane::arcane_aleph_hypre)
  set(FEMUTILS_HAS_PARALLEL_SOLVER TRUE)
  target_link_libraries(FemUtils PRIVATE Arcane::arcane_aleph_hypre)
  message(STATUS "HYPRE backend is available")
  set(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE TRUE)
  find_package(Hypre)
  # There is several target for Hypre depending on how Hypre is compiled
  # We try 'HYPRE::HYPRE' which is available with recents versions and
  # if it is not available, we try with arccon
  if (TARGET HYPRE::HYPRE)
    target_link_libraries(FemUtils PRIVATE HYPRE::HYPRE)
  elseif (TARGET arccon::Hypre)
    target_link_libraries(FemUtils PRIVATE arccon::Hypre)
  else()
    message(WARNING "Do not find a valid target for Hypre. This may produce link errors")
  endif()
  target_sources(FemUtils PRIVATE HypreDoFLinearSystem.cc)
endif()
endif()

if(NOT FEMUTILS_HAS_PARALLEL_SOLVER)
  message(STATUS "WARNING: No Arcane parallel solver is available. Parallel linear solver system is disabled")
endif()

# Add a variable in cache to check which parallel solvers are availables
set(FEMUTILS_HAS_PARALLEL_SOLVER ${FEMUTILS_HAS_PARALLEL_SOLVER} CACHE BOOL "Is Parallel solver available" FORCE)
set(FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS ${FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS} CACHE BOOL "Is Trilinos solver available" FORCE)
set(FEMUTILS_HAS_SOLVER_BACKEND_PETSC ${FEMUTILS_HAS_SOLVER_BACKEND_PETSC} CACHE BOOL "Is PETSc solver available" FORCE)
set(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE ${FEMUTILS_HAS_SOLVER_BACKEND_HYPRE} CACHE BOOL "Is Hypre solver available" FORCE)
