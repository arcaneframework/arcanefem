set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Poisson
  ${ACCELERATOR_SOURCES}
  ManufacturedSolutions.cc
  main.cc
  Fem_axl.h
  ManufacturedSolutions.h
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
    sphere_cut.hexa.msh

    circle_cut.msh
    circle_cut.quad.msh
    circle_cut.quad8.msh
    circle_cut.quad9.msh

    plancher.msh

    unit_square.msh
    unit_square.quad.msh

    mms.quad4_n4.msh
    mms.quad8_n4.msh
    mms.quad9_n4.msh

    3x3x3_cube_hexa8.msh
    3x3x3_cube_hexa20.msh
    3x3x3_cube_hexa27.msh
)
foreach (MESH_FILE IN LISTS MESH_FILES)
  file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach ()

enable_testing()
set(MODULE_NAME "poisson")

arcanefem_add_seq_test(NAME 2D_pntDirichlet
  COMMAND ./Poisson
  ARGS
    -A,//fem/solution-comparison-file=check/poisson_test_point_dirichlet_2D.txt
    inputs/perforatedSquare.pointDirichlet.2D.arc)

if (FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  arcanefem_add_seq_test(NAME 2D
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/petsc-flags=-log_view
      inputs/circle.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D.txt
      inputs/circle.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_quad
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad.txt
      inputs/circle.2D.quad.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D_quad.txt
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_seq_test(NAME 2D_scalar_neumann_quad
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_scalar_neumann_2D_quad.txt
      -A,//fem/boundary-conditions/neumann/value=2
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_seq_test(NAME 2D_quad8
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/circle_cut.quad8.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad8.txt
      inputs/circle.2D.quad.arc)
  arcanefem_add_seq_test(NAME 2D_quad9
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/circle_cut.quad9.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad9.txt
      inputs/circle.2D.quad.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad8
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad9
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad8
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad9
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_mms_quad4
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad4_n4.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME 2D_mms_quad8
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME 2D_mms_quad9
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME 3D_mms_hexa8
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa8.msh
      -A,//fem/manufactured-solution=quadratic
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME 3D_mms_hexa20
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa20.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME 3D_mms_hexa27
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa27.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME 3D_hexa20_cube
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa20.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa20.txt
      inputs/cube.3D.hexa.arc)
  arcanefem_add_seq_test(NAME 3D_hexa27_cube
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa27.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa27.txt
      inputs/cube.3D.hexa.arc)
  arcanefem_add_seq_test(NAME 3D
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      inputs/sphere.3D.arc)
  arcanefem_add_seq_test(NAME 3D_neumann
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_neumann_3D.txt
      inputs/sphere.neumann.3D.arc)
  arcanefem_add_seq_test(NAME 3D_hexa
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D_hexa.txt
      inputs/sphere.3D.hexa.arc)
  arcanefem_add_seq_test(NAME 3D_neumann_hexa
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_neumann_3D_hexa.txt
      inputs/sphere.neumann.3D.hexa.arc)
  arcanefem_add_seq_test(NAME 3D_scalar_neumann_hexa
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_scalar_neumann_3D_hexa.txt
      -A,//fem/boundary-conditions/neumann/value=3.1
      inputs/sphere.neumann.3D.hexa.arc)
  arcanefem_add_seq_test(NAME 2D_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/petsc-flags=-ksp_monitor
      inputs/circle.2D.arc)
  arcanefem_add_seq_test(NAME 2D_assembly_only
    COMMAND ./Poisson
    ARGS
      -A,//fem/solve-linear-system=false
      inputs/circle.2D.arc)
  arcanefem_add_seq_test(NAME 2D_parse_and_exit
    COMMAND ./Poisson
    ARGS
      -A,//fem/assemble-linear-system=false
      -A,//fem/solve-linear-system=false
      inputs/circle.2D.arc)
  arcanefem_add_seq_test(NAME manufactured_2D
    COMMAND ./Poisson
    ARGS
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME manufactured-poly_2D
    COMMAND ./Poisson
    ARGS
      -A,//fem/manufactured-solution=poly
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME manufactured-quadratic_2D
    COMMAND ./Poisson
    ARGS
      -A,//fem/manufactured-solution=quadratic
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME manufactured_2D_quad
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/unit_square.quad.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME manufactured-poly_2D_quad
    COMMAND ./Poisson
    ARGS
      -A,//fem/manufactured-solution=poly
      -A,//meshes/mesh/filename=meshes/unit_square.quad.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME manufactured-quadratic_2D_quad
    COMMAND ./Poisson
    ARGS
      -A,//fem/manufactured-solution=quadratic
      -A,//meshes/mesh/filename=meshes/unit_square.quad.msh
      inputs/manufactured.solution.arc)
  arcanefem_add_seq_test(NAME manufactured-quadratic_2D_quad8
    COMMAND ./Poisson
    ARGS
      -A,//fem/manufactured-solution=quadratic
      -A,//meshes/mesh/filename=meshes/circle_cut.quad8.msh
      inputs/manufactured.solution.arc)

  arcanefem_add_full_test(NAME 2D_bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=BSR
      inputs/circle.2D.arc)
  arcanefem_add_full_test(NAME 2D_af-bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=AF-BSR
      inputs/circle.2D.arc)
  arcanefem_add_full_test(NAME 2D_quad8_bsr
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/manufactured-solution=quadratic
      -A,//fem/manufactured-solution-tolerance=1e-10
      -A,//fem/matrix-format=BSR
      inputs/manufactured.solution.arc)
  arcanefem_add_full_test(NAME 2D_quad8_af-bsr
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/manufactured-solution=quadratic
      -A,//fem/manufactured-solution-tolerance=1e-10
      -A,//fem/matrix-format=AF-BSR
      inputs/manufactured.solution.arc)
  arcanefem_add_full_test(NAME 2D_quad9_bsr
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/manufactured-solution=quadratic
      -A,//fem/manufactured-solution-tolerance=1e-10
      -A,//fem/matrix-format=BSR
      inputs/manufactured.solution.arc)
  arcanefem_add_full_test(NAME 2D_quad9_af-bsr
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/manufactured-solution=quadratic
      -A,//fem/manufactured-solution-tolerance=1e-10
      -A,//fem/matrix-format=AF-BSR
      inputs/manufactured.solution.arc)
  arcanefem_add_full_test(NAME 3D_bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=BSR
      inputs/sphere.3D.arc)
  arcanefem_add_full_test(NAME 3D_af-bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=AF-BSR
      inputs/sphere.3D.arc)
  arcanefem_add_full_test(NAME 3D_hexa8_af-bsr
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa8.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa8.txt
      -A,//fem/matrix-format=AF-BSR
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_neumann_bsr-hexa
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_neumann_3D_hexa.txt
      -A,//fem/matrix-format=BSR
      inputs/sphere.neumann.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_scalar_neumann_bsr-hexa
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_scalar_neumann_3D_hexa.txt
      -A,//fem/matrix-format=BSR
      -A,//fem/boundary-conditions/neumann/value=3.1
      inputs/sphere.neumann.3D.hexa.arc)

  if (FEMUTILS_HAS_PARALLEL_SOLVER)
    arcanefem_add_mpi_test(NAME 2D
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
        inputs/circle.2D.arc)
    arcanefem_add_mpi_test(NAME 2D_quad
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad.txt
        inputs/circle.2D.quad.arc)
    arcanefem_add_mpi_test(NAME 2D_scalar_neumann_quad
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_scalar_neumann_2D_quad.txt
        -A,//fem/boundary-conditions/neumann/value=2
        inputs/circle.neumann.2D.quad.arc)
    arcanefem_add_mpi_test(NAME 3D
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
        inputs/sphere.3D.arc)
    arcanefem_add_mpi_test(NAME 3D_hexa
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D_hexa.txt
        inputs/sphere.3D.hexa.arc)
    arcanefem_add_mpi_test(NAME 3D_hexa8_cube
      COMMAND ./Poisson
      ARGS
        -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa8.msh
        -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa8.txt
        inputs/cube.3D.hexa.arc)
    arcanefem_add_mpi_test(NAME 3D_hexa20_cube
      COMMAND ./Poisson
      ARGS
        -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa20.msh
        -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa20.txt
        inputs/cube.3D.hexa.arc)
    arcanefem_add_mpi_test(NAME 3D_hexa27_cube
      COMMAND ./Poisson
      ARGS
        -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa27.msh
        -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa27.txt
        inputs/cube.3D.hexa.arc)
    arcanefem_add_mpi_test(NAME 3D_scalar_neumann_hexa
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_scalar_neumann_3D_hexa.txt
        -A,//fem/boundary-conditions/neumann/value=3.1
        inputs/sphere.neumann.3D.hexa.arc)
  endif ()

endif ()

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

  arcanefem_add_full_test(NAME 2D_pntDirichlet_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_point_dirichlet_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/perforatedSquare.pointDirichlet.2D.arc)
  arcanefem_add_full_test(NAME 2D_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.arc)
  arcanefem_add_full_test(NAME 2D_quad_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_neumann_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D_quad.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_scalar_neumann_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_scalar_neumann_2D_quad.txt
      -A,//fem/matrix-format=BSR
      -A,//fem/boundary-conditions/neumann/value=2
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad8_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad8.txt
      -A,//meshes/mesh/filename=meshes/circle_cut.quad8.msh
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad8_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad8.txt
      -A,//meshes/mesh/filename=meshes/circle_cut.quad8.msh
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad9_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/circle_cut.quad9.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad9.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad9_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/circle_cut.quad9.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad9.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad8_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad8_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad8_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad8_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=AF-BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad9_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad9_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad9_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad9_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=AF-BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_HYPRE_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_full_test(NAME 2D_quad_neumann_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D_quad.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_scalar_neumann_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_scalar_neumann_2D_quad.txt
      -A,//fem/matrix-format=AF-BSR
      -A,//fem/boundary-conditions/neumann/value=2
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_neumann_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.neumann.2D.arc)
  arcanefem_add_full_test(NAME 2D_bsr_atomicFree_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_2D}
      inputs/circle.2D.arc)
  arcanefem_add_full_test(NAME 3D_bsr_hypre_gpu
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/sphere.3D.arc)
  arcanefem_add_full_test(NAME 3D_neumann_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/sphere.neumann.3D.arc)
  arcanefem_add_full_test(NAME 3D_bsr_atomicFree_hypre_gpu
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/sphere.3D.arc)
  arcanefem_add_full_test(NAME 3D_hexa20_cube_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa20.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa20.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_hexa20_cube_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa20.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa20.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_hexa27_cube_bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa27.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa27.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_hexa27_cube_af-bsr_hypre
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa27.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa27.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_neumann_bsr-hexa_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_neumann_3D_hexa.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_3D}
      inputs/sphere.neumann.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_scalar_neumann_bsr-hexa_hypre
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_scalar_neumann_3D_hexa.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG_3D}
      -A,//fem/boundary-conditions/neumann/value=3.1
      inputs/sphere.neumann.3D.hexa.arc)
endif ()


if (FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  set(SOLVER_PETSC_CG_2D
    -A,//fem/linear-system/@name=PetscLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0)

  set(SOLVER_PETSC_CG_3D
    -A,//fem/linear-system/@name=PetscLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_full_test(NAME 2D_pntDirichlet_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_point_dirichlet_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/perforatedSquare.pointDirichlet.2D.arc)
  arcanefem_add_full_test(NAME 2D_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.arc)
  arcanefem_add_full_test(NAME 2D_quad_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_neumann_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D_quad.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_scalar_neumann_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_scalar_neumann_2D_quad.txt
      -A,//fem/matrix-format=BSR
      -A,//fem/boundary-conditions/neumann/value=2
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad8_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad8.txt
      -A,//meshes/mesh/filename=meshes/circle_cut.quad8.msh
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad8_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad8.txt
      -A,//meshes/mesh/filename=meshes/circle_cut.quad8.msh
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad9_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/circle_cut.quad9.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad9.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad9_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/circle_cut.quad9.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad9.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.quad.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad8_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad8_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad8_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad8_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad8_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad8.txt
      -A,//fem/matrix-format=AF-BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad9_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad9_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_neumann_quad9_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_seq_test(NAME 2D_vector_neumann_quad9_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/mms.quad9_n4.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_quadratic_neumann_2D_quad9.txt
      -A,//fem/matrix-format=AF-BSR
      "-A,//fem/boundary-conditions/neumann/value=2 0"
      ${SOLVER_PETSC_CG_2D}
      inputs/quadratic.neumann.2D.arc)
  arcanefem_add_full_test(NAME 2D_quad_neumann_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D_quad.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_quad_scalar_neumann_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_scalar_neumann_2D_quad.txt
      -A,//fem/matrix-format=AF-BSR
      -A,//fem/boundary-conditions/neumann/value=2
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.neumann.2D.quad.arc)
  arcanefem_add_full_test(NAME 2D_neumann_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.neumann.2D.arc)
  arcanefem_add_full_test(NAME 2D_bsr_atomicFree_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_2D}
      inputs/circle.2D.arc)
  arcanefem_add_full_test(NAME 3D_bsr_petsc_gpu
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/sphere.3D.arc)
  arcanefem_add_full_test(NAME 3D_neumann_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/sphere.neumann.3D.arc)
  arcanefem_add_full_test(NAME 3D_bsr_atomicFree_petsc_gpu
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/sphere.3D.arc)
  arcanefem_add_full_test(NAME 3D_hexa20_cube_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa20.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa20.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_hexa20_cube_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa20.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa20.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_hexa27_cube_bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa27.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa27.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_hexa27_cube_af-bsr_petsc
    COMMAND ./Poisson
    ARGS
      -A,//meshes/mesh/filename=meshes/3x3x3_cube_hexa27.msh
      -A,//fem/solution-comparison-file=check/poisson_test_ref_cube_3D_hexa27.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/cube.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_neumann_bsr-hexa_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_neumann_3D_hexa.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_3D}
      inputs/sphere.neumann.3D.hexa.arc)
  arcanefem_add_full_test(NAME 3D_scalar_neumann_bsr-hexa_petsc
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_scalar_neumann_3D_hexa.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_CG_3D}
      -A,//fem/boundary-conditions/neumann/value=3.1
      inputs/sphere.neumann.3D.hexa.arc)
endif ()

if (FEMUTILS_HAS_SOLVER_BACKEND_ALIEN)
  set(SOLVER_ALIEN_PETSC_CG
    -A,//fem/linear-system/@name=AlienLinearSystem
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0
    -A,//fem/linear-system/linear-solver/@name=PETScSolver
    -A,//fem/linear-system/linear-solver/solver/@name=CG
    -A,//fem/linear-system/linear-solver/solver/num-iterations-max=1000
    -A,//fem/linear-system/linear-solver/solver/stop-criteria-value=1e-15
    -A,//fem/linear-system/linear-solver/solver/preconditioner/@name=Jacobi
  )

  arcanefem_add_seq_test(NAME alien_petsc_2D_pntDirichlet_bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_point_dirichlet_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_ALIEN_PETSC_CG}
      inputs/perforatedSquare.pointDirichlet.2D.arc)

  arcanefem_add_seq_test(NAME alien_petsc_2D_bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_ALIEN_PETSC_CG}
      inputs/circle.2D.arc)

  arcanefem_add_seq_test(NAME alien_petsc_2D_neumann_bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_neumann_2D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_ALIEN_PETSC_CG}
      inputs/circle.neumann.2D.arc)

  arcanefem_add_seq_test(NAME alien_petsc_2D_bsr_atomicFree
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_ALIEN_PETSC_CG}
      inputs/circle.2D.arc)

  arcanefem_add_seq_test(NAME alien_petsc_3D_bsr_petsc_gpu
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_ALIEN_PETSC_CG}
      inputs/sphere.3D.arc)

  arcanefem_add_seq_test(NAME alien_petsc_3D_neumann_bsr
    COMMAND ./Poisson
    ARGS
      -A,//fem/matrix-format=BSR
      ${SOLVER_ALIEN_PETSC_CG}
      inputs/sphere.neumann.3D.arc)

  arcanefem_add_seq_test(NAME petsc_alien_3D_bsr_atomicFree_gpu
    COMMAND ./Poisson
    ARGS
      -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_ALIEN_PETSC_CG}
      inputs/sphere.3D.arc)

  if (FEMUTILS_HAS_PARALLEL_SOLVER)
    arcanefem_add_mpi_test(NAME alien_petsc_2D
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D.txt
        -A,//fem/matrix-format=BSR
        ${SOLVER_ALIEN_PETSC_CG}
        inputs/circle.2D.arc)

    arcanefem_add_mpi_test(NAME alien_petsc_2D_quad
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_2D_quad.txt
        -A,//fem/matrix-format=BSR
        ${SOLVER_ALIEN_PETSC_CG}
        inputs/circle.2D.quad.arc)

    arcanefem_add_mpi_test(NAME alien_petsc_2D_scalar_neumann_quad
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_circle_scalar_neumann_2D_quad.txt
        -A,//fem/boundary-conditions/neumann/value=2
        -A,//fem/matrix-format=BSR
        ${SOLVER_ALIEN_PETSC_CG}
        inputs/circle.neumann.2D.quad.arc)

    arcanefem_add_mpi_test(NAME alien_petsc_3D
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D.txt
        -A,//fem/matrix-format=BSR
        ${SOLVER_ALIEN_PETSC_CG}
        inputs/sphere.3D.arc)

    arcanefem_add_mpi_test(NAME alien_petsc_3D_hexa
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_3D_hexa.txt
        -A,//fem/matrix-format=BSR
        ${SOLVER_ALIEN_PETSC_CG}
        inputs/sphere.3D.hexa.arc)

    arcanefem_add_mpi_test(NAME alien_petsc_3D_scalar_neumann_hexa
      COMMAND ./Poisson
      ARGS
        -A,//fem/solution-comparison-file=check/poisson_test_ref_sphere_scalar_neumann_3D_hexa.txt
        -A,//fem/boundary-conditions/neumann/value=3.1
        -A,//fem/matrix-format=BSR
        ${SOLVER_ALIEN_PETSC_CG}
        inputs/sphere.neumann.3D.hexa.arc)
  endif ()
endif ()
