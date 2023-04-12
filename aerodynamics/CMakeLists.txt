cmake_minimum_required(VERSION 3.21)
project(aerodynamics LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(aerodynamics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(aerodynamics)
target_include_directories(aerodynamics PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(aerodynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Joukowski.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/NACA0012.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(aerodynamics PUBLIC FemUtils)

enable_testing()

add_test(NAME Joukowski COMMAND aerodynamics Test.Joukowski.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME Joukowski_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./aerodynamics Test.Joukowski.arc)
endif()
