set(ACCELERATOR_SOURCES
  FemModule.cc
  BlCsrBiliAssembly.cc
  CooGpuBiliAssembly.cc
  CooSortGpuBiliAssembly.cc
  CsrGpuBiliAssembly.cc
  NodeWiseCsrBiliAssembly.cc
)

add_executable(Testlab
  ${ACCELERATOR_SOURCES}
  FemModule.h
  LegacyBiliAssembly.cc
  CooBiliAssembly.cc
  CsrBiliAssembly.cc
  CooSortBiliAssembly.cc
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})

arcane_accelerator_add_to_target(Testlab)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Testlab)
target_include_directories(Testlab PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Testlab.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/L-shape.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/random.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/porous-medium.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/L-shape-3D.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/sphere_cut.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/circle_cut.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# For benchmarking
# Copy the benchmarking files in the binary directory
# The '/' after 'benchmarking' is needed because we want to copy the files
# inside the 'benchmarking' directory but not the directory itself.
file(COPY "benchmarking/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(Testlab PUBLIC FemUtils)

target_compile_definitions(Testlab PRIVATE $<$<BOOL:${ARCANE_HAS_ACCELERATOR}>:ARCANE_HAS_ACCELERATOR>)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

add_test(NAME [testlab]3D_Dirichlet COMMAND Testlab inputs/Test.L-shape.3D.arc)
if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [testlab]hypre_direct_2D_csr COMMAND Testlab inputs/Test.circle.2D.csr.arc)
  add_test(NAME [testlab]hypre_direct_2D_blcsr COMMAND Testlab inputs/Test.circle.2D.blcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_2D_blcsr_gpu COMMAND ./Testlab ARGS inputs/Test.circle.2D.blcsr.arc)
  add_test(NAME [testlab]hypre_direct_3D COMMAND Testlab inputs/Test.sphere.3D.arc)
  add_test(NAME [testlab]hypre_direct_3D_csr COMMAND Testlab inputs/Test.sphere.3D.csr.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [testlab]hypre_direct_2D_csr_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.circle.2D.csr.arc)
    add_test(NAME [testlab]hypre_direct_2D_blcsr_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.circle.2D.blcsr.arc)
    add_test(NAME [testlab]hypre_direct_3D_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.sphere.3D.arc)
    add_test(NAME [testlab]hypre_direct_3D_csr_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.sphere.3D.csr.arc)
  endif()
endif()


add_test(NAME [testlab]direct-solver COMMAND Testlab inputs/Test.direct-solver.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS)
  add_test(NAME [testlab]trilinos COMMAND Testlab inputs/Test.trilinos.arc)
endif()


if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [testlab]hypre COMMAND Testlab inputs/Test.hypre.arc)
  add_test(NAME [testlab]hypre_direct COMMAND Testlab inputs/Test.hypre_direct.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [testlab]hypre_direct_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.hypre_direct.arc)
    add_test(NAME [testlab]hypre_direct_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Testlab inputs/Test.hypre_direct.arc)
  endif()
endif()


if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [testlab]L-shape_2D COMMAND Testlab inputs/Test.L-shape.2D.arc)
  add_test(NAME [testlab]L-shape_2D_petsc_coo COMMAND Testlab inputs/Test.L-shape.2D.coo.arc)
  add_test(NAME [testlab]L-shape_2D_petsc_coo-sort COMMAND Testlab inputs/Test.L-shape.2D.coo-sort.arc)
  add_test(NAME [testlab]petsc COMMAND Testlab inputs/Test.petsc.arc)
  add_test(NAME [testlab]neumann COMMAND Testlab inputs/Test.neumann.arc)
  add_test(NAME [testlab]porous COMMAND Testlab inputs/Test.porous.arc)
  add_test(NAME [testlab]3D_Dirichlet_petsc_coo COMMAND Testlab inputs/Test.L-shape.3D.coo.arc)
  add_test(NAME [testlab]3D_Dirichlet_petsc_coo-sort COMMAND Testlab inputs/Test.L-shape.3D.coo-sort.arc)

  arcanefem_add_gpu_test(NAME [testlab]gpu COMMAND ./Testlab ARGS inputs/Test.petsc.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-sort-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-sort-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-sort-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-sort-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_nwcsr_gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.nwcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_nwcsr_gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.nwcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_csr-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.csr-gpu.arc)

  if (FEMUTILS_HAS_PARALLEL_SOLVER)
    add_test(NAME [testlab]L-shape_2D_petsc_coo_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.2D.coo.arc)
    add_test(NAME [testlab]L-shape_2D_petsc_coo-sort_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.2D.coo-sort.arc)
    add_test(NAME [testlab]3D_Dirichlet_petsc_coo_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.3D.coo.arc)
    add_test(NAME [testlab]3D_Dirichlet_petsc_coo-sort_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.3D.coo-sort.arc)

    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-sort-gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-sort-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-sort-gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-sort-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_nwcsr_gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.nwcsr.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_nwcsr_gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.nwcsr.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_csr-gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.csr-gpu.arc)
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_gpu COMMAND ./Testlab ARGS inputs/Test.hypre_direct.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D_gpu COMMAND ./Testlab ARGS inputs/Test.sphere.3D.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D_csr-gpu COMMAND ./Testlab ARGS inputs/Test.sphere.3D.csr-gpu.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER)
    arcanefem_add_gpu_test(NAME [testlab]hypre_direct_gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.hypre_direct.arc)
    arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D_gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.sphere.3D.arc)
    arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D_csr-gpu_2pe NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.sphere.3D.csr-gpu.arc)
  endif()
endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [testlab]L-shape_2D_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Testlab inputs/Test.L-shape.2D.arc)
endif()

add_test(NAME [testlab]2D_bsr COMMAND Testlab inputs/Test.L-shape.2D.bsr.arc)
arcanefem_add_gpu_test(NAME [testlab]2D_bsr_gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.bsr.arc)

add_test(NAME [testlab]3D_bsr COMMAND Testlab inputs/Test.sphere.3D.bsr.arc)
arcanefem_add_gpu_test(NAME [testlab]3D_bsr_gpu COMMAND ./Testlab ARGS inputs/Test.sphere.3D.bsr.arc)
