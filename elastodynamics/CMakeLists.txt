cmake_minimum_required(VERSION 3.21)
project(Elastodynamics LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(Elastodynamics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elastodynamics)
target_include_directories(Elastodynamics PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Elastodynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elastodynamics.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elastodynamics PUBLIC FemUtils)

enable_testing()
add_test(NAME elastodynamics COMMAND Elastodynamics Test.Elastodynamics.arc)
