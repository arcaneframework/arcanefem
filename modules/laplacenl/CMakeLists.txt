set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(LaplaceNL
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(LaplaceNL)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(LaplaceNL)

target_include_directories(LaplaceNL PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(LaplaceNL.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(LaplaceNL PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  ring.msh
  ring.quad.msh
  plancher.msh
  L-shape-3D.msh
  truncated_cube.msh
  truncated_cube.hexa.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/meshes/plancher.10k.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plancher.geo -format msh41 -o meshes/plancher.10k.msh -bin
    DEPENDS ${MSH_DIR}/plancher.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files_laplaceNL DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/meshes/plancher.10k.msh)
  add_dependencies(LaplaceNL gmsh_files_laplaceNL)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()

add_test(NAME [laplacenl]pointDirichlet COMMAND LaplaceNL
  -A,//fem/solution-comparison-file=check/test3_results.txt
  inputs/PointDirichlet.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [laplacenl]2D_Dirichlet COMMAND LaplaceNL
    inputs/ring.arc)
  add_test(NAME [laplacenl]2D_Dirichlet-Neumann COMMAND LaplaceNL
    inputs/ring.dirichlet-neumann.arc)
  add_test(NAME [laplacenl]2D_Dirichlet-Neumann_quad COMMAND LaplaceNL
    -A,//fem/solution-comparison-file=check/test_ring_quad.txt
    inputs/ring.dirichlet-neumann.quad.arc)
  add_test(NAME [laplacenl]3D_Dirichlet COMMAND LaplaceNL
    -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
    inputs/L-shape.3D.arc)
  add_test(NAME [laplacenl]3D_Dirichlet-Neumann COMMAND LaplaceNL
    inputs/L-shape.3D.dirichlet-neumann.arc)
  add_test(NAME [laplacenl]3D_Dirichlet-Neumann_hexa COMMAND LaplaceNL
    -A,//fem/solution-comparison-file=check/test_trucated-cube_hexa.txt
    inputs/truncated-cube.neumann.3D.hexa.arc)
  add_test(NAME [laplacenl]2D_Dirichlet_petsc_monitor COMMAND LaplaceNL
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/ring.arc)
  add_test(NAME [laplacenl]2D_parse_and_exit COMMAND LaplaceNL
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    inputs/ring.arc)

  arcanefem_add_gpu_test(NAME [laplacenl]2D_Dirichlet_bsr COMMAND ./LaplaceNL ARGS
    -A,//fem/matrix-format=BSR
    inputs/ring.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]2D_Dirichlet_af-bsr COMMAND ./LaplaceNL ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/ring.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_Dirichlet_bsr COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
    -A,//fem/matrix-format=BSR
    inputs/L-shape.3D.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_Dirichlet_af-bsr COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
    -A,//fem/matrix-format=AF-BSR
    inputs/L-shape.3D.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_pointBc_af-bsr COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test_3D_truncated-cube.txt
    -A,//fem/matrix-format=AF-BSR
    inputs/truncated-cube.3D.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_pointBc_neumann_af-bsr COMMAND ./LaplaceNL ARGS
    -A,//fem/matrix-format=AF-BSR
    -A,//fem/fp-tol=1e-3
    inputs/truncated-cube.neumann.3D.arc)

  set(SOLVER_PETSC_CG
    -A,//fem/linear-system/@name=PETScLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)
  arcanefem_add_gpu_test(NAME [laplacenl]pointDirichlet_bsr_petsc-direct COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test3_results.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_PETSC_CG}
    inputs/PointDirichlet.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]pointDirichlet_af-bsr_petsc-direct COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test3_results.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_PETSC_CG}
    inputs/PointDirichlet.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_Dirichlet_bsr_petsc-direct COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_PETSC_CG}
    inputs/L-shape.3D.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_Dirichlet_af-bsr_petsc-direct COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_PETSC_CG}
    inputs/L-shape.3D.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_pointBc_neumann_af-bsr_petsc-direct COMMAND ./LaplaceNL ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_PETSC_CG}
    inputs/truncated-cube.neumann.3D.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [laplacenl]2D_Dirichlet_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./LaplaceNL
      inputs/ring.arc)
    add_test(NAME [laplacenl]3D_Dirichlet_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./LaplaceNL
      -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
      inputs/L-shape.3D.arc)
    add_test(NAME [laplacenl]2D_Dirichlet-Neumann_quad_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./LaplaceNL
      inputs/ring.dirichlet-neumann.arc)
    add_test(NAME [laplacenl]3D_Dirichlet-Neumann_hexa_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./LaplaceNL
      -A,//fem/solution-comparison-file=check/test_trucated-cube_hexa.txt
      inputs/truncated-cube.neumann.3D.hexa.arc)

    if(FEMTEST_HAS_GMSH_TEST)
      add_test(NAME [laplacenl]pointDirichlet_10k_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./LaplaceNL
        -A,//meshes/mesh/filename=meshes/plancher.10k.msh
        -A,//fem/linear-system/@name=AlephLinearSystem
        -A,//fem/linear-system/solver-backend=petsc
        -A,//fem/petsc-flags=-ksp_monitor\ -ksp_rtol\ 10e-5
        -A,//fem/cross-validation=false
        inputs/PointDirichlet.arc)
    endif()
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)
  arcanefem_add_gpu_test(NAME [laplacenl]pointDirichlet_bsr_hypreDirect COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test3_results.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/PointDirichlet.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]pointDirichlet_af-bsr_hypreDirect COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test3_results.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/PointDirichlet.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_Dirichlet_bsr_hypreDirect COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/L-shape.3D.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_Dirichlet_af-bsr_hypreDirect COMMAND ./LaplaceNL ARGS
    -A,//fem/solution-comparison-file=check/test_3D_L-shape.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/L-shape.3D.arc)
  arcanefem_add_gpu_test(NAME [laplacenl]3D_pointBc_neumann_af-bsr_hypreDirect COMMAND ./LaplaceNL ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/truncated-cube.neumann.3D.arc)
endif()