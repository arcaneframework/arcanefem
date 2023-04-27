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
configure_file(Test.Bilaplacian.direct.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Bilaplacian.internal_pcg.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Bilaplacian.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bilap.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Bilaplacian PUBLIC FemUtils)

enable_testing()
add_test(NAME [bilaplacian]bilaplacian_direct COMMAND Bilaplacian Test.Bilaplacian.direct.arc)
add_test(NAME [bilaplacian]bilaplacian_internal_pcg COMMAND Bilaplacian Test.Bilaplacian.internal_pcg.arc)

#add_test(NAME [bilaplacian]bilaplacian COMMAND Bilaplacian Test.Bilaplacian.arc)
# If parallel part is available, add some tests
#if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
#  add_test(NAME [bilaplacian]bilaplacian_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Bilaplacian Test.Bilaplacian.arc)
#endif()
