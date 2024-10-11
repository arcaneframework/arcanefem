set(ACCELERATOR_SOURCES
  FemModule.cc
  CooBiliAssembly.cc
  CsrBiliAssembly.cc
  CooSortBiliAssembly.cc
  BlCsrBiliAssembly.cc
  CooGpuBiliAssembly.cc
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
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})

arcane_accelerator_add_to_target(Poisson)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Poisson)
target_include_directories(Poisson PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Poisson.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.2D.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.2D.coo.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.2D.coo-sort.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.3D.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.2D.coo-gpu.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.3D.coo-gpu.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.3D.coo.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.3D.coo-sort.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.sphere.3D.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.sphere.3D.csr.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.sphere.3D.csr.no-edge.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
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

target_compile_definitions(Poisson PRIVATE $<$<BOOL:${ARCANE_HAS_ACCELERATOR}>:ARCANE_HAS_ACCELERATOR>)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

add_test(NAME [poisson]3D_Dirichlet COMMAND Poisson Test.L-shape.3D.arc)
if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [poisson]hypre_direct_3D COMMAND Poisson Test.sphere.3D.arc)
  add_test(NAME [poisson]hypre_direct_3D_csr COMMAND Poisson Test.sphere.3D.csr.arc)
  add_test(NAME [poisson]hypre_direct_3D_csr_noedge COMMAND Poisson Test.sphere.3D.csr.no-edge.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [poisson]hypre_direct_3D_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Poisson Test.sphere.3D.arc)
    add_test(NAME [poisson]hypre_direct_3D_csr_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Poisson Test.sphere.3D.csr.arc)
    add_test(NAME [poisson]hypre_direct_3D_csr_noedge_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson Test.sphere.3D.csr.no-edge.arc)
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
  add_test(NAME [poisson]L-shape_2D_petsc_coo COMMAND Poisson Test.L-shape.2D.coo.arc)
  add_test(NAME [poisson]L-shape_2D_petsc_coo-sort COMMAND Poisson Test.L-shape.2D.coo-sort.arc)
  add_test(NAME [poisson]petsc COMMAND Poisson Test.petsc.arc)
  add_test(NAME [poisson]neumann COMMAND Poisson Test.neumann.arc)
  add_test(NAME [poisson]porous COMMAND Poisson Test.porous.arc)
  add_test(NAME [poisson]3D_Dirichlet_petsc_coo COMMAND Poisson Test.L-shape.3D.coo.arc)
  add_test(NAME [poisson]3D_Dirichlet_petsc_coo-sort COMMAND Poisson Test.L-shape.3D.coo-sort.arc)
endif()


arcanefem_add_gpu_test(NAME [poisson]gpu COMMAND ./Poisson ARGS Test.petsc.arc)
arcanefem_add_gpu_test(NAME [poisson]L-shape_2D_petsc_coo-gpu COMMAND ./Poisson ARGS Test.L-shape.2D.coo-gpu.arc)
arcanefem_add_gpu_test(NAME [poisson]L-shape_3D_petsc_coo-gpu COMMAND ./Poisson ARGS Test.L-shape.3D.coo-gpu.arc)
if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  arcanefem_add_gpu_test(NAME [poisson]hypre_direct_gpu COMMAND ./Poisson ARGS Test.hypre_direct.arc)
  arcanefem_add_gpu_test(NAME [poisson]hypre_direct_3D_gpu COMMAND ./Poisson ARGS Test.sphere.3D.arc)
  arcanefem_add_gpu_test(NAME [poisson]hypre_direct_3D_csr_noedge_gpu COMMAND ./Poisson ARGS Test.sphere.3D.csr.no-edge.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER)
    arcanefem_add_gpu_test(NAME [poisson]hypre_direct_gpu_2pe NB_MPI 2 COMMAND ./Poisson ARGS Test.hypre_direct.arc)
    arcanefem_add_gpu_test(NAME [poisson]hypre_direct_3D_gpu_2pe NB_MPI 2 COMMAND ./Poisson ARGS Test.sphere.3D.arc)
  endif()
endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [poisson]L-shape_2D_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson Test.L-shape.2D.arc)
endif()
