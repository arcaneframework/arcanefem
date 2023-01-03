cmake_minimum_required(VERSION 3.21)
project(Bilaplacian LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(Bilaplacian
  Fem1Module.cc
  main.cc
  Fem1_axl.h
)

arcane_generate_axl(Fem1)
arcane_add_arcane_libraries_to_target(Bilaplacian)
target_include_directories(Bilaplacian PUBLIC . ../fem1 ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FemTest1.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(FemTest1.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(bilap.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Bilaplacian PUBLIC FemUtils)

enable_testing()
add_test(NAME bilaplacian COMMAND Bilaplacian FemTest1.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME bilaplacian_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Bilaplacian FemTest1.arc)
endif()
