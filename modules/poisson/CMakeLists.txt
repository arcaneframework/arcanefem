set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Poisson
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Poisson)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Poisson)
target_include_directories(Poisson PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Poisson.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Poisson PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
    sphere_cut.msh
    circle_cut.msh
    plancher.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

add_test(NAME [poisson]2D_pntDirichlet COMMAND Poisson
    inputs/perforatedSquare.pointDirichlet.2D.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [poisson]2D COMMAND Poisson
    inputs/circle.2D.arc)
  add_test(NAME [poisson]2D_neumann COMMAND Poisson
    inputs/circle.neumann.2D.arc)
  add_test(NAME [poisson]3D COMMAND Poisson
    inputs/sphere.3D.arc)
  add_test(NAME [poisson]3D_neumann COMMAND Poisson
    inputs/sphere.neumann.3D.arc)
  add_test(NAME [poisson]2D_petsc COMMAND Poisson
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/circle.2D.arc)
  add_test(NAME [poisson]2D_assembly_only COMMAND Poisson
    -A,//fem/solve-linear-system=false
    -A,//fem/cross-validation=false
    inputs/circle.2D.arc)
  add_test(NAME [poisson]2D_no_cross_validation COMMAND Poisson
    -A,//fem/cross-validation=false
    inputs/circle.2D.arc)
  add_test(NAME [poisson]2D_parse_and_exit COMMAND Poisson
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    -A,//fem/cross-validation=false
    inputs/circle.2D.arc)

  arcanefem_add_gpu_test(NAME [poisson]2D_bsr COMMAND Poisson ARGS
    -A,//fem/matrix-format=BSR
    inputs/circle.2D.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_af-bsr COMMAND Poisson ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/circle.2D.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_bsr COMMAND Poisson ARGS
    -A,//fem/matrix-format=BSR
    inputs/sphere.3D.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_af-bsr COMMAND Poisson ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/sphere.3D.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER)
    add_test(NAME [poisson]2D_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Poisson
      inputs/circle.2D.arc)
    add_test(NAME [poisson]3D_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Poisson
      inputs/sphere.3D.arc)
  endif()

endif()

if (FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG_2D
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0)

  set(SOLVER_HYPRE_CG_3D
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/preconditioner=amg
    -A,//fem/linear-system/amg-threshold=0.55
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_gpu_test(NAME [poisson]2D_pntDirichlet_bsr_hypre COMMAND Poisson ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/perforatedSquare.pointDirichlet.2D.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_bsr_hypre COMMAND Poisson ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/circle.2D.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_neumann_bsr_hypre COMMAND Poisson ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/circle.neumann.2D.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_bsr_atomicFree_hypre COMMAND Poisson ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/circle.2D.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_bsr_hypre_gpu COMMAND Poisson ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_3D}
    inputs/sphere.3D.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_neumann_bsr_hypre COMMAND Poisson ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_3D}
    inputs/sphere.neumann.3D.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_bsr_atomicFree_hypre_gpu COMMAND Poisson ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG_3D}
    inputs/sphere.3D.arc)
endif()