cmake_minimum_required(VERSION 3.21)
project(Bilaplacian LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(Bilaplacian
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Bilaplacian)
target_include_directories(Bilaplacian PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Bilaplacian.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.internal_pcg.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.direct.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bilap.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Bilaplacian PUBLIC FemUtils)

enable_testing()

add_test(NAME [bilaplacian]direct_solver COMMAND Bilaplacian Test.direct.arc)
add_test(NAME [bilaplacian]internal_PCG_solver COMMAND Bilaplacian Test.internal_pcg.arc)
