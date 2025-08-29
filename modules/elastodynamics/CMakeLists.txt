set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Elastodynamics
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elastodynamics)
target_include_directories(Elastodynamics PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Elastodynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elastodynamics PUBLIC FemUtils)

# Copy the data files in the binary directory
file(COPY "data" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  bar_dynamic.msh
  bar_dynamic_quad.msh
  bar_dynamic_3D.msh
  bar_dynamic_3Dhexa.msh
  semi-circle.msh
  truncated_cube.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [elastodynamics]2D_bar COMMAND Elastodynamics
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 2\ -ksp_view\ -ksp_monitor
    inputs/bar.arc)

  add_test(NAME [elastodynamics]2D_bar_RowElimination COMMAND Elastodynamics
    -A,//fem/boundary-conditions/dirichlet[1]/enforce-Dirichlet-method=RowElimination
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 2\ -ksp_view\ -ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-10
    inputs/bar.arc)

  add_test(NAME [elastodynamics]2D_bar_RowColumnElimination COMMAND Elastodynamics
    -A,//fem/boundary-conditions/dirichlet[1]/enforce-Dirichlet-method=RowColumnElimination
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 2\ -ksp_view\ -ksp_monitor\  -ksp_type\ gmres\ -pc_type\ bjacobi
    inputs/bar.arc)

  add_test(NAME [elastodynamics]2D_bar_quad COMMAND Elastodynamics
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 2\ -ksp_view\ -ksp_monitor
    inputs/bar.quad.arc)

  add_test(NAME [elastodynamics]3D_bar COMMAND Elastodynamics
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_view\ -ksp_monitor
    inputs/bar.3D.arc)

  add_test(NAME [elastodynamics]3D_bar_RowColumnElimination_bjacobi COMMAND Elastodynamics
    -A,//fem/boundary-conditions/dirichlet[1]/enforce-Dirichlet-method=RowColumnElimination
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-10
    inputs/bar.3D.arc)

  add_test(NAME [elastodynamics]3D_bar_hexa COMMAND Elastodynamics
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_view\ -ksp_monitor
    inputs/bar.3D.hexa.arc)

  add_test(NAME [elastodynamics]2D_transient_traction COMMAND Elastodynamics
    inputs/bar.transient-traction.arc)

  add_test(NAME [elastodynamics]Dirichlet_pointBc_Penalty COMMAND Elastodynamics
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 2\ -ksp_view\ -ksp_monitor\ -ksp_type\ cg\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15
    inputs/semi-circle.pointBC.arc)

  add_test(NAME [elastodynamics]2D_pointBc_RowElimination COMMAND Elastodynamics
    -A,//fem/boundary-conditions/dirichlet[1]/enforce-Dirichlet-method=RowElimination
    -A,//fem/boundary-conditions/dirichlet-point[1]/enforce-Dirichlet-method=RowElimination
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 2\ -ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-10
    inputs/semi-circle.pointBC.arc)

  add_test(NAME [elastodynamics]3D_truncated-cube_pointBc_Penalty COMMAND Elastodynamics
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_view\ -ksp_monitor\ -ksp_type\ cg\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15
    inputs/truncated-cube.pointBC.arc)

  add_test(NAME [elastodynamics]3D_truncated-cube_pointBc_RowColElim COMMAND Elastodynamics
    -A,//fem/boundary-conditions/dirichlet[1]/enforce-Dirichlet-method=RowColumnElimination
    -A,//fem/boundary-conditions/dirichlet-point[1]/enforce-Dirichlet-method=RowColumnElimination
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_view\ -ksp_monitor\ -ksp_type\ cg\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15
    inputs/truncated-cube.pointBC.arc)

  add_test(NAME [elastodynamics]3D_truncated-cube_pointBc_RowElim COMMAND Elastodynamics
    -A,//fem/boundary-conditions/dirichlet[1]/enforce-Dirichlet-method=RowElimination
    -A,//fem/boundary-conditions/dirichlet-point[1]/enforce-Dirichlet-method=RowElimination
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_view\ -ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15
    inputs/truncated-cube.pointBC.arc)

  add_test(NAME [elastodynamics]3D_transient_traction COMMAND Elastodynamics
    inputs/bar.3D.transient-traction.arc)

  add_test(NAME [elastodynamics]2D_constant_traction_and_damping COMMAND Elastodynamics
    inputs/bar.damping.arc)

  add_test(NAME [elastodynamics]2D_time-discretization_Galpha COMMAND Elastodynamics
    inputs/bar.Galpha.arc)

  add_test(NAME [elastodynamics]parse_and_exit COMMAND Elastodynamics
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    -A,//fem/cross-validation=false
    inputs/bar.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]2D_bar_petsc_bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=BSR
    -A,//fem/petsc-flags=-ksp_view\ -ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-10
    inputs/bar.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]3D_bar_petsc_bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=BSR
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_view\ -ksp_monitor
    inputs/bar.3D.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]2D_bar_petsc_af-bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=AF-BSR
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 2\ -ksp_view\ -ksp_monitor
    inputs/bar.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]3D_bar_petsc_af-bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=AF-BSR
    -A,//fem/petsc-flags=-mat_type\ baij\ -mat_block_size\ 3\ -ksp_view\ -ksp_monitor
    inputs/bar.3D.arc)

endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
  -A,//fem/linear-system/@name=HypreLinearSystem
  -A,//fem/linear-system/solver=cg
  -A,//fem/linear-system/atol=1e-9
  -A,//fem/linear-system/rtol=0)

  add_test(NAME [elastodynamics]Dirichlet_traction_bodyforce COMMAND Elastodynamics
    inputs/bar.dirichlet.traction.bodyforce.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]2D_bar_hypre_bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]3D_bar_hypre_bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.3D.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]2D_bar_bsr_af-bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.arc)

  arcanefem_add_gpu_test(NAME [elastodynamics]3D_bar_hypre_af-bsr COMMAND  ./Elastodynamics ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/bar.3D.arc)

  add_test(NAME [elastodynamics]3D_truncated-cube_hypre_pointBc_Penalty COMMAND Elastodynamics
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/truncated-cube.pointBC.arc)
endif()
