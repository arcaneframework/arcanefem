cmake_minimum_required(VERSION 3.21)
project(Poisson LANGUAGES C CXX)

find_package(Arcane REQUIRED)
find_library(CUSPARSE_LIB cusparse HINTS ${CUDA_TOOLKIT_ROOT_DIR}/lib64)

option(REGISTER_TIME "Enable time registration" ON)
option(USE_COO_GPU "Use the CSR ArcaneFEM format instead of the linear hashmap, but computing on GPU" ON)
option(USE_COO "Use the CSR ArcaneFEM format instead of the linear hashmap" ON)
option(USE_CUSPARSE "Use the cusparse library instead of the linear hashmap" OFF)
option(USE_LEGACY "Use the linear hashmap (original code)" ON)
arcane_accelerator_enable()
arcane_accelerator_add_source_files(
  FemModule.cc
)


add_executable(Poisson
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_to_target(Poisson)

target_link_libraries(Poisson PRIVATE ${CUSPARSE_LIB})

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Poisson)
target_include_directories(Poisson PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Poisson.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.poisson.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.poisson.direct.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.poisson.neumann.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.poisson.trilinos.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.poisson.hypre.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.poisson.petsc.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/L-shape.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/random.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Poisson PUBLIC FemUtils)

target_compile_definitions(Poisson PRIVATE $<$<BOOL:${REGISTER_TIME}>:REGISTER_TIME>)
target_compile_definitions(Poisson PRIVATE $<$<BOOL:${USE_CUSPARSE}>:USE_CUSPARSE>)
target_compile_definitions(Poisson PRIVATE $<$<BOOL:${USE_COO}>:USE_COO>)
target_compile_definitions(Poisson PRIVATE $<$<BOOL:${USE_COO_GPU}>:USE_COO_GPU>)
target_compile_definitions(Poisson PRIVATE $<$<BOOL:${USE_LEGACY}>:USE_LEGACY>)

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/L-100-shape.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/L-shape.geo -setnumber rfactor 100 -format msh41 -L-100-shape.msh
    DEPENDS ${MSH_DIR}/L-shape.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
add_custom_target(gmsh_files_poisson DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/L-100-shape.msh)
  add_dependencies(Poisson gmsh_files_poisson)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

add_test(NAME [poisson]poisson COMMAND Poisson Test.poisson.arc)
add_test(NAME [poisson]poisson_direct COMMAND Poisson Test.poisson.direct.arc)
add_test(NAME [poisson]poisson_neumann COMMAND Poisson Test.poisson.neumann.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS)
  add_test(NAME [poisson]poisson_trilinos COMMAND Poisson Test.poisson.trilinos.arc)
endif()
if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [poisson]poisson_hypre COMMAND Poisson Test.poisson.hypre.arc)
endif()
if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [poisson]poisson_petsc COMMAND Poisson Test.poisson.petsc.arc)
endif()


# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [poisson]poisson_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson Test.poisson.arc)
endif()
