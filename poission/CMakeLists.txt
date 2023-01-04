cmake_minimum_required(VERSION 3.21)
project(FemTest LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(FemTest
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(FemTest)
target_include_directories(FemTest PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FemTest.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.heterogeneous.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.heterogeneous.10k.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.10k.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.quad4.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.laplace.PointDirichlet.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(multi-material.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(plancher.quad4.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(FemTest PUBLIC FemUtils)

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
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/multi-material.10k.msh
    COMMAND ${GMSH} -2 ${CMAKE_CURRENT_SOURCE_DIR}/multi-mat.geo -format msh41 -o multi-material.10k.msh
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/multi-material.10k.msh ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh)
  add_dependencies(FemTest gmsh_files)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()

add_test(NAME conduction COMMAND FemTest Test.conduction.arc)
add_test(NAME conduction_heterogeneous COMMAND FemTest Test.conduction.heterogeneous.arc)
add_test(NAME conduction_quad COMMAND FemTest Test.conduction.quad4.arc)
add_test(NAME laplace_pointDirichlet COMMAND FemTest Test.laplace.PointDirichlet.arc)
if(FEMTEST_HAS_GMSH_TEST)
  add_test(NAME conduction_10k COMMAND ./FemTest Test.conduction.10k.arc)
  add_test(NAME conduction_heterogeneous_10k COMMAND ./FemTest Test.conduction.heterogeneous.10k.arc)
endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME conduction_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./FemTest Test.conduction.arc)
  add_test(NAME conduction_heterogeneous_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./FemTest Test.conduction.heterogeneous.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME conduction_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./FemTest Test.conduction.10k.arc)
    add_test(NAME conduction_heterogeneous_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./FemTest Test.conduction.heterogeneous.10k.arc)
  endif()
endif()
