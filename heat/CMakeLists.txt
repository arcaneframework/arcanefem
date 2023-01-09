cmake_minimum_required(VERSION 3.21)
project(conduction LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(conduction
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(conduction)
target_include_directories(conduction PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FemTest.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(plate.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(conduction PUBLIC FemUtils)

enable_testing()

add_test(NAME conduction COMMAND conduction Test.conduction.arc)
