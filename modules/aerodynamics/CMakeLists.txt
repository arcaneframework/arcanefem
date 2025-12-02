set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(aerodynamics
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(aerodynamics)
target_include_directories(aerodynamics PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(aerodynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(aerodynamics PUBLIC FemUtils)

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the check folder containing the arc files
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  NACA0012.msh
  NACA0012.quad.msh
  aerodynamics_3d_coarse.msh
  aerodynamics_3d_coarse.hexa.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [aerodynamics]2d_Joukowski COMMAND aerodynamics
    -A,//fem/solution-comparison-file=check/test_2d.txt
    inputs/Joukowski.arc)
  add_test(NAME [aerodynamics]2d_Joukowski_quad COMMAND aerodynamics
    -A,//fem/solution-comparison-file=check/test_2d.quad.txt
    inputs/Joukowski.quad.arc)
  add_test(NAME [aerodynamics]parse_and_exit COMMAND aerodynamics
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    inputs/Joukowski.arc)

  arcanefem_add_gpu_test(NAME [aerodynamics]2d_Joukowski_bsr COMMAND ./aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_2d.txt
    -A,//fem/matrix-format=BSR
    inputs/Joukowski.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]2d_Joukowski_af-bsr COMMAND ./aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_2d.txt
    -A,//fem/matrix-format=AF-BSR
    inputs/Joukowski.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]3d_Joukowski_bsr COMMAND ./aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_3d.txt
    -A,//fem/matrix-format=BSR
    inputs/Joukowski_3d.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]3d_Joukowski_af-bsr COMMAND ./aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_3d.txt
    -A,//fem/matrix-format=AF-BSR
    inputs/Joukowski_3d.arc)

  set(SOLVER_PETSC_CG_2D
    -A,//fem/linear-system/@name=PETScLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0)

  set(SOLVER_PETSC_CG_3D
    -A,//fem/linear-system/@name=PETScLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-15
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_gpu_test(NAME [aerodynamics]2d_Joukowski_bsr_petsc COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_2d.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_PETSC_CG_2D}
    inputs/Joukowski.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]2d_Joukowski_af_petsc COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_2d.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_PETSC_CG_2D}
    inputs/Joukowski.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]3d_Joukowski_bsr_petsc COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_3d.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_PETSC_CG_3D}
    inputs/Joukowski_3d.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]3d_Joukowski_bsr_petsc_gpu COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_3d.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_PETSC_CG_3D}
    inputs/Joukowski_3d.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [aerodynamics]2d_Joukowski_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./aerodynamics
      -A,//fem/solution-comparison-file=check/test_2d.txt
      inputs/Joukowski.arc)
    add_test(NAME [aerodynamics]2d_Joukowski_quad_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./aerodynamics
      -A,//fem/solution-comparison-file=check/test_2d.quad.txt
      inputs/Joukowski.quad.arc)
    add_test(NAME [aerodynamics]3d_Joukowski_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./aerodynamics
      -A,//fem/solution-comparison-file=check/test_3d.txt
      inputs/Joukowski_3d.arc)
    add_test(NAME [aerodynamics]3d_Joukowski_hexa_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./aerodynamics
      -A,//fem/solution-comparison-file=check/test_3d.hexa.txt
      inputs/Joukowski_3d.hexa.arc)
    add_test(NAME [aerodynamics]3d_Joukowski_bjacobi_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./aerodynamics
      -A,//fem/solution-comparison-file=check/test_3d.txt
      -A,//fem/petsc-flags=-ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15
      inputs/Joukowski_3d.arc)
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

  arcanefem_add_gpu_test(NAME [aerodynamics]2d_Joukowski_bsr_hypre COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_2d.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/Joukowski.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]2d_Joukowski_af_hypre COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_2d.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG_2D}
    inputs/Joukowski.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]3d_Joukowski_bsr_hypre COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_3d.txt
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG_3D}
    inputs/Joukowski_3d.arc)
  arcanefem_add_gpu_test(NAME [aerodynamics]3d_Joukowski_bsr_hypre_gpu COMMAND aerodynamics ARGS
    -A,//fem/solution-comparison-file=check/test_3d.txt
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG_3D}
    inputs/Joukowski_3d.arc)
endif()
