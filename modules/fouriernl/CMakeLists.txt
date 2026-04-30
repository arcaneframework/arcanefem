set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(FourierNL
  FemModule.h
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(FourierNL)
target_include_directories(FourierNL PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FourierNL.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(FourierNL PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
    unit_square.msh
    unit_square.coarse.msh
    unit_square.quad.msh
    unit_cube.msh
    unit_cube.hexa-non-cartesian.msh)

foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [fouriernl]2d_nonlinear_conduction COMMAND FourierNL
          -A,//fem/solution-comparison-file=check/conduction_tria.txt
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.arc)
  add_test(NAME [fouriernl]2d_nonlinear_conduction_qdot_zero COMMAND FourierNL
          -A,//fem/qdot=0.0
          -A,//fem/solution-comparison-file=check/conduction_tria.txt
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.arc)
  add_test(NAME [fouriernl]2d_nonlinear_conduction_quad COMMAND FourierNL
          -A,//fem/solution-comparison-file=check/conduction_quad.txt
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.quad.arc)
  add_test(NAME [fouriernl]2d_nonlinear_conduction_coarse COMMAND FourierNL
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.coarse.arc)
  add_test(NAME [fouriernl]3d_nonlinear_conduction_hexa COMMAND FourierNL
          -A,//fem/solution-comparison-file=check/conduction_hexa.txt
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.3d.nonlinear.conduction.hexa.arc)
  add_test(NAME [fouriernl]3d_nonlinear_conduction_hexa_cartesian COMMAND FourierNL
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.3d.nonlinear.conduction.hexa.cartesian.arc)
  add_test(NAME [fouriernl]3d_nonlinear_conduction COMMAND FourierNL
          -A,//fem/solution-comparison-file=check/conduction_tetra.txt
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.3d.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]2d_conduction_petsc_bsr COMMAND ./FourierNL ARGS
          -A,//fem/solution-comparison-file=check/conduction_tria.txt
          -A,//fem/matrix-format=BSR
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]2d_conduction_petsc_af-bsr COMMAND ./FourierNL ARGS
          -A,//fem/solution-comparison-file=check/conduction_tria.txt
          -A,//fem/matrix-format=AF-BSR
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]3d_conduction_petsc_bsr COMMAND ./FourierNL ARGS
          -A,//fem/matrix-format=BSR
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.3d.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]3d_conduction_petsc_af-bsr COMMAND ./FourierNL ARGS
          -A,//fem/matrix-format=AF-BSR
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.3d.nonlinear.conduction.arc)

  set(SOLVER_PETSC_CG
    -A,//fem/linear-system/@name=PETScLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_gpu_test(NAME [fouriernl]2d_conduction_petsc-direct_bsr COMMAND ./FourierNL ARGS
          -A,//fem/matrix-format=BSR
          ${SOLVER_PETSC_CG}
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]2d_conduction_petsc-direct_af-bsr COMMAND ./FourierNL ARGS
          -A,//fem/matrix-format=AF-BSR
          ${SOLVER_PETSC_CG}
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]3d_conduction_petsc-direct_bsr COMMAND ./FourierNL ARGS
          -A,//fem/matrix-format=BSR
          ${SOLVER_PETSC_CG}
          -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
          inputs/Test.3d.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]3d_conduction_petsc-direct_af-bsr COMMAND ./FourierNL ARGS
          -A,//fem/matrix-format=AF-BSR
          ${SOLVER_PETSC_CG}
          inputs/Test.3d.nonlinear.conduction.arc)


  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [fouriernl]2d_nonlinear_conduction_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./FourierNL
            -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
            inputs/Test.nonlinear.conduction.arc)
    add_test(NAME [fouriernl]2d_nonlinear_conduction_quad_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./FourierNL
            -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
            inputs/Test.nonlinear.conduction.quad.arc)
    add_test(NAME [fouriernl]2d_nonlinear_conduction_coarse_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./FourierNL
            -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
            inputs/Test.nonlinear.conduction.coarse.arc)
    add_test(NAME [fouriernl]3d_nonlinear_conduction_hexa_cartesian_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./FourierNL
            -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
            inputs/Test.3d.nonlinear.conduction.hexa.cartesian.arc
            -A,//meshes/mesh/generator/nb-part-x=2)
    add_test(NAME [fouriernl]3d_nonlinear_conduction_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./FourierNL
            -A,//fem/petsc-flags=-ksp_initial_guess_nonzero
            inputs/Test.3d.nonlinear.conduction.arc)
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_gpu_test(NAME [fouriernl]2d_conduction_hypre_bsr COMMAND ./FourierNL ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/Test.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]2d_conduction_hypre_af-bsr COMMAND ./FourierNL ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/Test.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]3d_conduction_hypre_bsr COMMAND ./FourierNL ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/Test.3d.nonlinear.conduction.arc)

  arcanefem_add_gpu_test(NAME [fouriernl]3d_conduction_hypre_af-bsr COMMAND ./FourierNL ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/Test.3d.nonlinear.conduction.arc)

endif()