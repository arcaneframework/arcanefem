set (USE_CUSPARSE_ADD OFF)

if (ARCANE_HAS_CUDA)
  if (TARGET CUDA::cusparse)
    message(STATUS "cuSPARSE is available. Enabling cusparse assembly")
    set(USE_CUSPARSE_ADD ON)
  else()
    message(STATUS "cuSPARSE is not available")
  endif()
endif()

set(ACCELERATOR_SOURCES
  FemModule.cc
  CooBiliAssembly.cc
  CsrBiliAssembly.cc
  CooSortBiliAssembly.cc
  BlCsrBiliAssembly.cc
  CsrGpuBiliAssembly.cc
  NodeWiseCsrBiliAssembly.cc
)

add_executable(Poisson
  ${ACCELERATOR_SOURCES}
  FemModule.h
  LegacyBiliAssembly.cc
  FemModule.cc
  main.cc
  Fem_axl.h
  CusparseBiliAssembly.hxx
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})

arcane_accelerator_add_to_target(Poisson)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Poisson)
target_include_directories(Poisson PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Poisson.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.2D.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.3D.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.sphere.3D.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.sphere.3D.csr.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.direct-solver.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.neumann.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.porous.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.trilinos.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.hypre.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.hypre_direct.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.petsc.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/L-shape.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/random.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/porous-medium.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/L-shape-3D.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/sphere_cut.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)


target_link_libraries(Poisson PUBLIC FemUtils)

if (USE_CUSPARSE_ADD)
  target_compile_definitions(Poisson PRIVATE USE_CUSPARSE_ADD)
  target_link_libraries(Poisson PRIVATE CUDA::cusparse)
endif()
target_compile_definitions(Poisson PRIVATE $<$<BOOL:${ARCANE_HAS_ACCELERATOR}>:ARCANE_HAS_ACCELERATOR>)

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/L-100-shape.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/L-shape.geo -setnumber rfactor 100 -format msh41 -o ${CMAKE_CURRENT_BINARY_DIR}/L-100-shape.msh
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

add_test(NAME [poisson]3D_Dirichlet COMMAND Poisson Test.L-shape.3D.arc)
if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [poisson]hypre_direct_3D COMMAND Poisson Test.sphere.3D.arc)
  add_test(NAME [poisson]hypre_direct_3D_csr COMMAND Poisson Test.sphere.3D.csr.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [poisson]hypre_direct_3D_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Poisson Test.sphere.3D.arc)
    add_test(NAME [poisson]hypre_direct_3D_csr_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Poisson Test.sphere.3D.csr.arc)
  endif()
endif()


add_test(NAME [poisson]direct-solver COMMAND Poisson Test.direct-solver.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS)
  add_test(NAME [poisson]trilinos COMMAND Poisson Test.trilinos.arc)
endif()


if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [poisson]hypre COMMAND Poisson Test.hypre.arc)
  add_test(NAME [poisson]hypre_direct COMMAND Poisson Test.hypre_direct.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [poisson]hypre_direct_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Poisson Test.hypre_direct.arc)
    add_test(NAME [poisson]hypre_direct_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson Test.hypre_direct.arc)
  endif()
endif()


if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [poisson]L-shape_2D COMMAND Poisson Test.L-shape.2D.arc)
  add_test(NAME [poisson]petsc COMMAND Poisson Test.petsc.arc)
  add_test(NAME [poisson]neumann COMMAND Poisson Test.neumann.arc)
  add_test(NAME [poisson]porous COMMAND Poisson Test.porous.arc)
endif()


arcanefem_add_gpu_test(NAME [poisson]gpu COMMAND ./Poisson ARGS Test.petsc.arc)
if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  arcanefem_add_gpu_test(NAME [poisson]hypre_direct_gpu COMMAND ./Poisson ARGS Test.hypre_direct.arc)
  arcanefem_add_gpu_test(NAME [poisson]hypre_direct_3D_gpu COMMAND ./Poisson ARGS Test.sphere.3D.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER)
    arcanefem_add_gpu_test(NAME [poisson]hypre_direct_gpu_2pe NB_MPI 2 COMMAND ./Poisson ARGS Test.hypre_direct.arc)
    arcanefem_add_gpu_test(NAME [poisson]hypre_direct_3D_gpu_2pe NB_MPI 2 COMMAND ./Poisson ARGS Test.sphere.3D.arc)
  endif()
endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [poisson]L-shape_2D_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson Test.L-shape.2D.arc)
endif()
