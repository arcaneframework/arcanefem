set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(heat
  ${ACCELERATOR_SOURCES}
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(heat)
target_include_directories(heat PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Heat.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  plate.msh
  plate.quad.msh
  multi-material.msh
  multi-material.quad.msh
  truncated_cube.msh
  truncated_cube.hexa.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

target_link_libraries(heat PUBLIC FemUtils)

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/meshes/plate.fine.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plate.fine.geo -format msh41 -o meshes/plate.fine.msh -bin
    DEPENDS ${MSH_DIR}/plate.fine.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files_heat DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/meshes/plate.fine.msh)
  add_dependencies(heat gmsh_files_heat)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()


if(FEMUTILS_HAS_SOLVER_BACKEND_ALIEN)

  # <linear-system name="AlienLinearSystem">
  #   <linear-solver name="PETScSolver">
  #     <solver name="CG">
  #       <num-iterations-max>1000</num-iterations-max>
  #       <stop-criteria-value>1e-15</stop-criteria-value>
  #       <preconditioner name="Jacobi"></preconditionner>
  #     </solver>
  #   </linear-solver>
  # </linear-system>

  set(SOLVER_ALIEN_PETSC
          -A,//fem/linear-system/@name=AlienLinearSystem
          -A,//fem/linear-system/linear-solver/@name=PETScSolver
          -A,//fem/linear-system/linear-solver/solver/@name=CG
          -A,//fem/linear-system/linear-solver/solver/num-iterations-max=1000
          -A,//fem/linear-system/linear-solver/solver/stop-criteria-value=1e-15
          -A,//fem/linear-system/linear-solver/solver/preconditioner/@name=Jacobi
  )

  add_test(NAME [heat]alien_petsc_conduction COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.arc)

  add_test(NAME [heat]alien_petsc_conduction_quad COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction.quad.txt
          -A,//fem/matrix-format=AF-BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.quad.arc)

  add_test(NAME [heat]alien_petsc_conduction_RowElimination_Dirichlet COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction.txt
          -A,//fem/boundary-conditions/dirichlet/enforce-Dirichlet-method=RowElimination
          -A,//fem/linear-system/solver-backend=petsc
          -A,//fem/linear-system/solver-method=gmres
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.arc)

  add_test(NAME [heat]alien_petsc_conduction_quad_RowElimination_Dirichlet COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction.quad.txt
          -A,//fem/boundary-conditions/dirichlet/enforce-Dirichlet-method=RowElimination
          -A,//fem/linear-system/solver-backend=petsc
          -A,//fem/linear-system/solver-method=gmres
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.quad.arc)

  add_test(NAME [heat]alien_petsc_conduction_multi-material COMMAND heat
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.multi-material.arc)

  add_test(NAME [heat]alien_petsc_conduction_multi-material_quad COMMAND heat
          inputs/conduction.multi-material.quad.arc)

  add_test(NAME [heat]alien_petsc_conduction_convection COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction_convection.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.convection.arc)

  add_test(NAME [heat]alien_petsc_conduction_convection_quad COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction_convection.quad.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.convection.quad.arc)

  add_test(NAME [heat]alien_petsc_conduction_neumann_pointBC COMMAND heat
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.neumann.point-dirichlet.arc)

  add_test(NAME [heat]alien_petsc_conduction_neumann_pointBC_quad COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction_neumann_pointBC.quad.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.neumann.point-dirichlet.quad.arc)

  add_test(NAME [heat]alien_petsc_conduction_petsc_monitor COMMAND heat
          -A,//fem/solution-comparison-file=check/2d_conduction.txt
          -A,//fem/petsc-flags=-ksp_monitor
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.arc)

  add_test(NAME [heat]alien_petsc_2d_parse_and_exit COMMAND heat
          -A,//fem/assemble-linear-system=false
          -A,//fem/solve-linear-system=false
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/conduction.arc)

  add_test(NAME [heat]alien_petsc_3d_conduction COMMAND heat
          -A,//fem/solution-comparison-file=check/3d_conduction.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/3d_conduction.arc)

  add_test(NAME [heat]alien_petsc_3d_conduction_hexa COMMAND heat
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/3d_conduction.hexa.arc)

  add_test(NAME [heat]alien_petsc_3d_conduction_neumann COMMAND heat
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/3d_conduction.neumann.arc)

  add_test(NAME [heat]alien_petsc_3d_conduction_neumann_hexa COMMAND heat
          -A,//fem/solution-comparison-file=check/3d_conduction_neumann.hexa.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/3d_conduction.neumann.hexa.arc)

  add_test(NAME [heat]alien_petsc_3d_conduction_convection_hexa COMMAND heat
          -A,//fem/solution-comparison-file=check/3d_conduction_convection.hexa.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/3d_conduction.convection.hexa.arc)

  add_test(NAME [heat]alien_petsc_3d_conduction_convection_pointBC COMMAND heat
          -A,//fem/solution-comparison-file=check/3d_conduction_convection_pointBC.txt
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/3d_conduction.pointBc.convection.arc)

  add_test(NAME [heat]alien_petsc_3d_conduction_convection_pointBC_rowElim COMMAND heat
          -A,//fem/solution-comparison-file=check/3d_conduction_convection_pointBC.txt
          -A,//fem/boundary-conditions/dirichlet/enforce-Dirichlet-method=RowElimination
          -A,//fem/boundary-conditions/dirichlet-point/enforce-Dirichlet-method=RowElimination
          -A,//fem/linear-system/solver-backend=petsc
          -A,//fem/linear-system/solver-method=gmres
          -A,//fem/matrix-format=BSR
          ${SOLVER_ALIEN_PETSC}
          inputs/3d_conduction.pointBc.convection.arc)

endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [heat]conduction COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction.txt
    inputs/conduction.arc)

  add_test(NAME [heat]conduction_quad COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction.quad.txt
    inputs/conduction.quad.arc)

  add_test(NAME [heat]conduction_RowElimination_Dirichlet COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction.txt
    -A,//fem/boundary-conditions/dirichlet/enforce-Dirichlet-method=RowElimination
    -A,//fem/linear-system/solver-backend=petsc
    -A,//fem/linear-system/solver-method=gmres
    inputs/conduction.arc)

  add_test(NAME [heat]conduction_quad_RowElimination_Dirichlet COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction.quad.txt
    -A,//fem/boundary-conditions/dirichlet/enforce-Dirichlet-method=RowElimination
    -A,//fem/linear-system/solver-backend=petsc
    -A,//fem/linear-system/solver-method=gmres
    inputs/conduction.quad.arc)

  add_test(NAME [heat]conduction_multi-material COMMAND heat
    inputs/conduction.multi-material.arc)

  add_test(NAME [heat]conduction_multi-material_quad COMMAND heat
    inputs/conduction.multi-material.quad.arc)

  add_test(NAME [heat]conduction_convection COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction_convection.txt
    inputs/conduction.convection.arc)

  add_test(NAME [heat]conduction_convection_quad COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction_convection.quad.txt
    inputs/conduction.convection.quad.arc)

  add_test(NAME [heat]conduction_neumann_pointBC COMMAND heat
    inputs/conduction.neumann.point-dirichlet.arc)

  add_test(NAME [heat]conduction_neumann_pointBC_quad COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction_neumann_pointBC.quad.txt
    inputs/conduction.neumann.point-dirichlet.quad.arc)

  add_test(NAME [heat]conduction_petsc_monitor COMMAND heat
    -A,//fem/solution-comparison-file=check/2d_conduction.txt
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/conduction.arc)

  add_test(NAME [heat]2d_parse_and_exit COMMAND heat
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    inputs/conduction.arc)

  add_test(NAME [heat]3d_conduction COMMAND heat
    -A,//fem/solution-comparison-file=check/3d_conduction.txt
    inputs/3d_conduction.arc)

  add_test(NAME [heat]3d_conduction_hexa COMMAND heat
    inputs/3d_conduction.hexa.arc)

  add_test(NAME [heat]3d_conduction_neumann COMMAND heat
    inputs/3d_conduction.neumann.arc)

  add_test(NAME [heat]3d_conduction_neumann_hexa COMMAND heat
    -A,//fem/solution-comparison-file=check/3d_conduction_neumann.hexa.txt
    inputs/3d_conduction.neumann.hexa.arc)

  add_test(NAME [heat]3d_conduction_convection_hexa COMMAND heat
    -A,//fem/solution-comparison-file=check/3d_conduction_convection.hexa.txt
    inputs/3d_conduction.convection.hexa.arc)

  add_test(NAME [heat]3d_conduction_convection_pointBC COMMAND heat
    -A,//fem/solution-comparison-file=check/3d_conduction_convection_pointBC.txt
    inputs/3d_conduction.pointBc.convection.arc)

  add_test(NAME [heat]3d_conduction_convection_pointBC_rowElim COMMAND heat
    -A,//fem/solution-comparison-file=check/3d_conduction_convection_pointBC.txt
    -A,//fem/boundary-conditions/dirichlet/enforce-Dirichlet-method=RowElimination
    -A,//fem/boundary-conditions/dirichlet-point/enforce-Dirichlet-method=RowElimination
    -A,//fem/linear-system/solver-backend=petsc
    -A,//fem/linear-system/solver-method=gmres
    inputs/3d_conduction.pointBc.convection.arc)

  arcanefem_add_gpu_test(NAME [heat]conduction_petsc_bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/2d_conduction.txt
    -A,//fem/matrix-format=BSR
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [heat]conduction_petsc_af-bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/2d_conduction.txt
    -A,//fem/matrix-format=AF-BSR
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [heat]3d_conduction_petsc_bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/3d_conduction.txt
    -A,//fem/matrix-format=BSR
    inputs/3d_conduction.arc)

  arcanefem_add_gpu_test(NAME [heat]3d_conduction_petsc_af-bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/3d_conduction.txt
    -A,//fem/matrix-format=AF-BSR
    inputs/3d_conduction.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [heat]conduction_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/2d_conduction.txt
      inputs/conduction.arc)

    add_test(NAME [heat]conduction_quad_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/2d_conduction.quad.txt
      inputs/conduction.quad.arc)

    add_test(NAME [heat]3d_conduction_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/3d_conduction.txt
      inputs/3d_conduction.arc)

    add_test(NAME [heat]3d_conduction_hexa_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      inputs/3d_conduction.hexa.arc)

    add_test(NAME [heat]3d_conduction__convection_hexa_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/3d_conduction_convection.hexa.txt
      inputs/3d_conduction.convection.hexa.arc)

    add_test(NAME [heat]conduction_neumann_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      inputs/conduction.neumann.point-dirichlet.arc)

    add_test(NAME [heat]conduction_neumann_quad_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/2d_conduction_neumann_pointBC.quad.txt
      inputs/conduction.neumann.point-dirichlet.quad.arc)

    add_test(NAME [heat]conduction_RowElimination_Dirichlet_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/2d_conduction.txt
      -A,//fem/boundary-conditions/dirichlet/enforce-Dirichlet-method=RowElimination
      -A,//fem/linear-system/solver-backend=petsc
      -A,//fem/linear-system/solver-method=gmres
      inputs/conduction.arc)

    add_test(NAME [heat]conduction_convection_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/2d_conduction_convection.txt
      inputs/conduction.convection.arc)

    add_test(NAME [heat]conduction_convection_quad_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      -A,//fem/solution-comparison-file=check/2d_conduction_convection.quad.txt
      inputs/conduction.convection.quad.arc)

    add_test(NAME [heat]conduction_neumann_pointBC_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat
      inputs/conduction.neumann.point-dirichlet.arc)

    if(FEMTEST_HAS_GMSH_TEST)
      add_test(NAME [heat]conduction_fine_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.fine.arc)
      add_test(NAME [heat]conduction_convection_fine_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.convection.fine.arc)
    endif()
  endif()

endif()

if (FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG_2D
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_gpu_test(NAME [heat]conduction_hypre_bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/2d_conduction.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [heat]conduction_hypre_af-bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/2d_conduction.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [heat]3d_conduction_hypre_bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/3d_conduction.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/3d_conduction.arc)

  arcanefem_add_gpu_test(NAME [heat]3d_conduction_hypre_af-bsr COMMAND heat ARGS
    -A,//fem/solution-comparison-file=check/3d_conduction.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/3d_conduction.arc)
endif()