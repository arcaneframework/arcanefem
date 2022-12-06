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
)

arcane_generate_axl(Fem1)
arcane_add_arcane_libraries_to_target(FemTest1)
target_include_directories(FemTest1 PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FemTest1.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(FemTest1.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
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

enable_testing()
add_test(NAME femtest1 COMMAND FemTest1 FemTest1.arc)

# If parallel part is available, add some tests
if(FEMTEST_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  set(MPIEXEC_EXECUTABLE /ccc/spack/links4/gcc9-mpich3/bin/mpiexec)

  add_test(NAME femtest1_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./FemTest1 FemTest1.arc)
endif()
