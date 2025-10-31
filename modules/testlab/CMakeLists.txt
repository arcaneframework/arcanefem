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

add_test(NAME [testlab]direct-solver COMMAND Testlab inputs/Test.direct-solver.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_TRILINOS)
  add_test(NAME [testlab]trilinos COMMAND Testlab inputs/Test.trilinos.arc)
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [testlab]circle_2D_petsc_traction_dok COMMAND Testlab inputs/Test.circle.2D.trac.dok.arc)
  add_test(NAME [testlab]circle_2D_petsc_x-traction_dok COMMAND Testlab inputs/Test.circle.2D.trac-x.dok.arc)
  add_test(NAME [testlab]circle_2D_petsc_y-traction_dok COMMAND Testlab inputs/Test.circle.2D.trac-y.dok.arc)
  add_test(NAME [testlab]circle_2D_petsc_vect-traction_dok COMMAND Testlab inputs/Test.circle.2D.trac2.dok.arc)
  add_test(NAME [testlab]L-shape_2D_petsc_coo COMMAND Testlab inputs/Test.L-shape.2D.coo.arc)
  add_test(NAME [testlab]L-shape_2D_petsc_coo-sort COMMAND Testlab inputs/Test.L-shape.2D.coo-sort.arc)
  add_test(NAME [testlab]neumann COMMAND Testlab inputs/Test.neumann.arc)
  add_test(NAME [testlab]porous COMMAND Testlab inputs/Test.porous.arc)
  add_test(NAME [testlab]3D_Dirichlet COMMAND Testlab inputs/Test.L-shape.3D.arc)
  add_test(NAME [testlab]3D_Dirichlet_petsc_coo COMMAND Testlab inputs/Test.L-shape.3D.coo.arc)
  add_test(NAME [testlab]3D_Dirichlet_petsc_coo-sort COMMAND Testlab inputs/Test.L-shape.3D.coo-sort.arc)
  add_test(NAME [testlab]2D_all-dok-matrix-format_cache-warm COMMAND  Testlab -A,DOK=TRUE -A,cache_warming=2 inputs/Test.petsc.arc)
  add_test(NAME [testlab]3D_all-dok-matrix-format_cache-warm COMMAND  Testlab -A,DOK=TRUE -A,cache_warming=2 inputs/Test.L-shape.3D.arc)
  add_test(NAME [testlab]2D_all-coo-matrix-formats_cache-warm COMMAND  Testlab -A,COO=TRUE,S-COO=TRUE,COO_GPU=TRUE,S-COO_GPU=TRUE -A,cache_warming=2 inputs/Test.porous.arc)
  add_test(NAME [testlab]3D_all-coo-matrix-formats_cache-warm COMMAND  Testlab -A,COO=TRUE,S-COO=TRUE,COO_GPU=TRUE,S-COO_GPU=TRUE -A,cache_warming=2 inputs/Test.L-shape.3D.arc)
  add_test(NAME [testlab]2D_all-csr-matrix-formats_cache-warm COMMAND  Testlab -A,CSR=TRUE,CSR_GPU=TRUE,AF-CSR_GPU=TRUE,BL-CSR_GPU=TRUE -A,cache_warming=2 inputs/Test.porous.arc)
  add_test(NAME [testlab]3D_all-csr-matrix-formats_cache-warm COMMAND  Testlab -A,CSR=TRUE,CSR_GPU=TRUE,AF-CSR_GPU=TRUE,BL-CSR_GPU=TRUE -A,cache_warming=2 inputs/Test.L-shape.3D.arc)
  add_test(NAME [testlab]2D_all-bsr-matrix-formats_cache-warm COMMAND  Testlab -A,BSR_GPU=TRUE,AF-BSR_GPU=TRUE -A,cache_warming=2 inputs/Test.L-shape.2D.bsr.arc)
  add_test(NAME [testlab]3D_all-bsr-matrix-formats_cache-warm COMMAND  Testlab -A,BSR_GPU=TRUE,AF-BSR_GPU=TRUE -A,cache_warming=2 inputs/Test.sphere.3D.bsr.arc)

  arcanefem_add_gpu_test(NAME [testlab]petsc_direct COMMAND  Testlab ARGS inputs/Test.petsc_direct.arc)
  arcanefem_add_gpu_test(NAME [testlab]petsc_direct_3D COMMAND  Testlab ARGS inputs/Test.petsc_direct-3D.arc)

  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc COMMAND ./Testlab ARGS inputs/Test.petsc.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-sort-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-sort-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_nwcsr COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.nwcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_csr-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.csr-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_bsr COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_af-bsr COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.atomic-free.arc)
  arcanefem_add_gpu_test(NAME [testlab]sphere_3D_petsc_bsr COMMAND Testlab ARGS inputs/Test.sphere.3D.bsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]sphere_3D_petsc_af-bsr COMMAND Testlab ARGS inputs/Test.sphere.3D.bsr.atomic-free.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-sort-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-sort-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_nwcsr COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.nwcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-gpu COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-gpu.arc)

  if (FEMUTILS_HAS_PARALLEL_SOLVER)
    add_test(NAME [testlab]circle_2D_petsc_vect-traction_dok_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.circle.2D.trac2.dok.arc)
    add_test(NAME [testlab]circle_2D_petsc_vect-traction-x_dok_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.circle.2D.trac-x.dok.arc)
    add_test(NAME [testlab]circle_2D_petsc_vect-traction-y_dok_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.circle.2D.trac-y.dok.arc)
    add_test(NAME [testlab]L-shape_2D_petsc_coo_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.2D.coo.arc)
    add_test(NAME [testlab]L-shape_2D_petsc_coo-sort_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.2D.coo-sort.arc)
    add_test(NAME [testlab]3D_Dirichlet_petsc_coo_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.3D.coo.arc)
    add_test(NAME [testlab]3D_Dirichlet_petsc_coo-sort_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.L-shape.3D.coo-sort.arc)

    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-gpu_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-gpu_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_coo-sort-gpu_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.coo-sort-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_petsc_coo-sort-gpu_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.coo-sort-gpu.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_nwcsr_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.nwcsr.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_3D_nwcsr_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.3D.nwcsr.arc)
    arcanefem_add_gpu_test(NAME [testlab]L-shape_2D_petsc_csr-gpu_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.L-shape.2D.csr-gpu.arc)
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [testlab]L-shape_2D_hypre COMMAND Testlab inputs/Test.hypre.arc)
  add_test(NAME [testlab]hypre_direct_2D_csr COMMAND Testlab inputs/Test.circle.2D.csr.arc)
  add_test(NAME [testlab]hypre_direct_trac_2D_csr COMMAND Testlab inputs/Test.circle.2D.trac.csr.arc)
  add_test(NAME [testlab]hypre_direct_vect-trac_2D_csr COMMAND Testlab inputs/Test.circle.2D.trac2.csr.arc)
  add_test(NAME [testlab]hypre_direct_3D_csr COMMAND Testlab inputs/Test.sphere.3D.csr.arc)

  arcanefem_add_gpu_test(NAME [testlab]hypre_direct COMMAND ./Testlab ARGS inputs/Test.hypre_direct.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D COMMAND ./Testlab ARGS inputs/Test.sphere.3D.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D_csr-gpu COMMAND ./Testlab ARGS inputs/Test.sphere.3D.csr-gpu.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct_gmres COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.gmres.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct_fgmres COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.fgmres.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct_bicgstab COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.bicgstab.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct_cg_bjac COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.cg.bjacobi.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct_gmres_bjac COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.gmres.bjacobi.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct_fgmres_bjac COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.fgmres.bjacobi.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_bsr_hypre_direct_bicgstab_bjac COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.hypre.bicgstab.bjacobi.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_2D_blcsr COMMAND ./Testlab ARGS inputs/Test.circle.2D.blcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_vect-trac_2D_blcsr COMMAND ./Testlab ARGS inputs/Test.circle.2D.trac2.blcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_trac-x_2D_blcsr COMMAND ./Testlab ARGS inputs/Test.circle.2D.trac-x.blcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]hypre_direct_trac-y_2D_blcsr COMMAND ./Testlab ARGS inputs/Test.circle.2D.trac-y.blcsr.arc)
  arcanefem_add_gpu_test(NAME [testlab]2D_af-bsr_hypre COMMAND Testlab ARGS inputs/Test.L-shape.2D.bsr.atomic-free.hypre.arc)
  arcanefem_add_gpu_test(NAME [testlab]3D_bsr_hypre COMMAND Testlab ARGS inputs/Test.sphere.3D.bsr.hypre.arc)
  arcanefem_add_gpu_test(NAME [testlab]3D_af-bsr_hypre COMMAND Testlab ARGS inputs/Test.sphere.3D.bsr.atomic-free.hypre.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER)
    add_test(NAME [testlab]hypre_direct_2D_blcsr_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.circle.2D.blcsr.arc)
    add_test(NAME [testlab]hypre_direct_3D_csr_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.sphere.3D.csr.arc)
    add_test(NAME [testlab]hypre_direct_2D_csr_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Testlab inputs/Test.circle.2D.csr.arc)
    add_test(NAME [testlab]L-2D_bsr_hypre_direct_gmres_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Testlab inputs/Test.L-shape.2D.bsr.hypre.gmres.arc)

    arcanefem_add_gpu_test(NAME [testlab]hypre_direct_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.hypre_direct.arc)
    arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.sphere.3D.arc)
    arcanefem_add_gpu_test(NAME [testlab]hypre_direct_3D_csr-gpu_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.sphere.3D.csr-gpu.arc)

    arcanefem_add_gpu_test(NAME [testlab]petsc_direct_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.petsc_direct.arc)
    arcanefem_add_gpu_test(NAME [testlab]petsc_direct_4p NB_MPI 4 COMMAND ./Testlab ARGS inputs/Test.petsc_direct.arc)
    arcanefem_add_gpu_test(NAME [testlab]petsc_direct_3D_2p NB_MPI 2 COMMAND ./Testlab ARGS inputs/Test.petsc_direct-3D.arc)
    arcanefem_add_gpu_test(NAME [testlab]petsc_direct_3D_4p NB_MPI 4 COMMAND ./Testlab ARGS inputs/Test.petsc_direct-3D.arc)
  endif()
endif()