set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Elasticity
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Elasticity)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elasticity)

target_include_directories(Elasticity PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})

configure_file(Elasticity.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elasticity PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
    bar.msh
    bar_dynamic_3D.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  set(SOLVER_PETSC_GMRES
    -A,//fem/petsc-flags=-ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15)

  add_test(NAME [elasticity]Dirichlet COMMAND Elasticity
    inputs/bar.2D.arc)

  add_test(NAME [elasticity]Dirichlet_traction COMMAND Elasticity
    inputs/bar.2D.traction.arc)

  add_test(NAME [elasticity]Dirichlet_pointBC COMMAND Elasticity
    inputs/bar.2D.PointDirichlet.arc)

  add_test(NAME [elasticity]Dirichlet_pointBC_via_RowElimination COMMAND Elasticity
    -A,//fem/enforce-Dirichlet-method=RowElimination
    ${SOLVER_PETSC_GMRES}
    inputs/bar.2D.PointDirichlet.arc)

  add_test(NAME [elasticity]Dirichlet_pointBC_via_RowColElimination COMMAND Elasticity
    -A,//fem/enforce-Dirichlet-method=RowColumnElimination
    inputs/bar.2D.PointDirichlet.arc)

  add_test(NAME [elasticity]Dirichlet_via_RowElimination COMMAND Elasticity
    -A,//fem/enforce-Dirichlet-method=RowElimination
    ${SOLVER_PETSC_GMRES}
    inputs/bar.2D.arc)

  add_test(NAME [elasticity]Dirichlet_via_RowColElimination COMMAND Elasticity
    -A,//fem/enforce-Dirichlet-method=RowColumnElimination
    inputs/bar.2D.arc)

  add_test(NAME [elasticity]Dirichlet_traction_bodyforce COMMAND Elasticity
    inputs/bar.2D.traction.bodyforce.arc)

  add_test(NAME [elasticity]bar_parse_and_exit COMMAND Elasticity
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    -A,//fem/cross-validation=false
    inputs/bar.2D.arc)

  add_test(NAME [elasticity]3D_Dirichlet_bodyforce_dok COMMAND Elasticity
    -A,//fem/matrix-format=DOK
    inputs/bar.2D.PointDirichlet.arc)

  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_pointBC_bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    inputs/bar.3D.Dirichlet.bodyForce.arc)

  arcanefem_add_gpu_test(NAME [elasticity]3D_Dirichlet_bodyforce_bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    inputs/bar.3D.Dirichlet.bodyForce.arc)

  arcanefem_add_gpu_test(NAME [elasticity]3D_Dirichlet_bodyforce_af-bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/bar.3D.Dirichlet.bodyForce.arc)

  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    inputs/bar.2D.traction.arc)

  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_bodyforce_bsr COMMAND ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    inputs/bar.2D.traction.bodyforce.arc)

  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_bodyforce_af-bsr COMMAND ./Elasticity ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/bar.2D.traction.bodyforce.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER)
    add_test(NAME [elasticity]Dirichlet_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity
      inputs/bar.2D.arc)
    add_test(NAME [elasticity]Dirichlet_via_RowElimination_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity
      -A,//fem/enforce-Dirichlet-method=RowElimination
      ${SOLVER_PETSC_GMRES}
      inputs/bar.2D.arc)
    add_test(NAME [elasticity]Dirichlet_via_RowColElimination_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity
      -A,//fem/enforce-Dirichlet-method=RowColumnElimination
      inputs/bar.2D.arc)
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)

    set(SOLVER_HYPRE_BICGSTAB
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=bicgstab)

  arcanefem_add_gpu_test(NAME [elasticity]3D_Dirichlet_bodyforce_traction_hypre_bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.3D.Dirichlet.traction.bodyForce.arc)
  arcanefem_add_gpu_test(NAME [elasticity]3D_Dirichlet_bodyforce_traction_hypre_af-bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.3D.Dirichlet.traction.bodyForce.arc)
  arcanefem_add_gpu_test(NAME [elasticity]3D_Dirichlet_bodyforce_hypre_bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.3D.Dirichlet.bodyForce.arc)
  arcanefem_add_gpu_test(NAME [elasticity]3D_Dirichlet_bodyforce_hypre_af-bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.3D.Dirichlet.bodyForce.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_pointBC_RowElimination_bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_BICGSTAB}
    -A,//fem/enforce-Dirichlet-method=RowElimination
    inputs/bar.2D.PointDirichlet.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_pointBC_RowElimination_af-bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_BICGSTAB}
    -A,//fem/enforce-Dirichlet-method=RowElimination
    inputs/bar.2D.PointDirichlet.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_pointBC_hypre_bsr COMMAND  ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    -A,//fem/enforce-Dirichlet-method=Penalty
    inputs/bar.2D.PointDirichlet.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_bsr_hypre COMMAND ./Elasticity ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.2D.traction.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_bodyforce_af-bsr_hypre COMMAND ./Elasticity ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.2D.traction.bodyforce.arc)
endif()