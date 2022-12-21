cmake_minimum_required(VERSION 3.21)
project(FemTest1 LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(FemTest1
  FemUtils.h
  FemUtils.cc
  Fem1Module.cc
  main.cc
  Fem1_axl.h
  FemLinearSystem.h
  FemLinearSystem.cc
  AlephFemLinearSystem.cc
  ../fem1/FemLinearSystem2.h
  ../fem1/FemLinearSystem2.cc
  ../fem1/AlephFemLinearSystem2.cc
  ../fem1/FemDoFsOnNodes.h
  ../fem1/FemDoFsOnNodes.cc
)

arcane_generate_axl(Fem1)
arcane_add_arcane_libraries_to_target(FemTest1)
target_include_directories(FemTest1 PUBLIC . ../fem1 ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FemTest1.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(FemTest1.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(FemTest1.10k.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(FemTest1 PRIVATE Arcane::arcane_aleph)
if (TARGET Arcane::arcane_aleph_hypre)
  set(FEMTEST_HAS_PARALLEL_SOLVER TRUE)
  target_link_libraries(FemTest1 PRIVATE Arcane::arcane_aleph_hypre)
else()
  message(STATUS "WARNING: Target Arcane::arcane_aleph_hypre is not available. Parallel linear solver system is disabled")
endif()

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh
    COMMAND ${GMSH} -2 ${CMAKE_CURRENT_SOURCE_DIR}/plancher.geo -format msh41 -o plancher.10k.msh
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh)
  add_dependencies(FemTest1 gmsh_files)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()
add_test(NAME femtest1 COMMAND FemTest1 FemTest1.arc)

# If parallel part is available, add some tests
if(FEMTEST_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME femtest1_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./FemTest1 FemTest1.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME femtest1_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./FemTest1 FemTest1.10k.arc)
  endif()
endif()
