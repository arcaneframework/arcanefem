cmake_minimum_required(VERSION 3.21)
project(Laplace LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(Laplace
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Laplace)
target_include_directories(Laplace PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Laplace.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.laplace.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.laplace.PointDirichlet.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.laplace.PointDirichlet.10K.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/ring.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Laplace PUBLIC FemUtils)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plancher.geo -format msh41 -o plancher.10k.msh
    DEPENDS ${MSH_DIR}/plancher.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files_laplace DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh)
  add_dependencies(Laplace gmsh_files_laplace)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()

add_test(NAME [laplace]laplace COMMAND Laplace Test.laplace.arc)
add_test(NAME [laplace]laplace_pointDirichlet COMMAND Laplace Test.laplace.PointDirichlet.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [laplace]laplace_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace Test.laplace.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME [laplace]laplace_pointDirichlet_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace Test.laplace.PointDirichlet.10K.arc)
  endif()
endif()
