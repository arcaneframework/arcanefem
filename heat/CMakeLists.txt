cmake_minimum_required(VERSION 3.21)
project(heat LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(heat
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(heat)
target_include_directories(heat PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Heat.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.convection.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(plate.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(heat PUBLIC FemUtils)

enable_testing()

add_test(NAME conduction COMMAND heat Test.conduction.arc)
add_test(NAME conduction_convection COMMAND heat Test.conduction.convection.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME conduction_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.arc)
  add_test(NAME conduction_convection_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.convection.arc)
endif()
