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
)

arcane_generate_axl(Fem1)
arcane_add_arcane_libraries_to_target(FemTest1)
target_include_directories(FemTest1 PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FemTest1.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(FemTest1.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
