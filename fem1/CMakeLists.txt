cmake_minimum_required(VERSION 3.21)
project(FemTest1 LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(FemTest1 Fem1Module.cc main.cc Fem1_axl.h)

arcane_generate_axl(Fem1)
arcane_add_arcane_libraries_to_target(FemTest1)
target_include_directories(FemTest1 PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FemTest1.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
