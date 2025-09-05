set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Soildynamics
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Soildynamics)
target_include_directories(Soildynamics PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Soildynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Soildynamics PUBLIC FemUtils)

# Copy the data files in the binary directory
file(COPY "data"  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  bar_dynamic.msh
  semi-circle-soil.msh
  square_double-couple.msh
  soil_2d.msh
  cube_double_couple_3d.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [soildynamics]2d_semi-circle_const_traction COMMAND Soildynamics
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/constant-traction.arc)
  add_test(NAME [soildynamics]2d_semi-circle_const_traction_pointbc COMMAND Soildynamics
    -A,//fem/petsc-flags=-ksp_view
    inputs/constant-traction.pointbc.arc)
  add_test(NAME [soildynamics]2d_semi-circle_const-traction_pointbc_rowelim COMMAND Soildynamics
    -A,//fem/petsc-flags=-ksp_view
    -A,//fem/boundary-conditions/dirichlet-point[1]/enforce-Dirichlet-method=RowElimination
    inputs/constant-traction.pointbc.arc)
  add_test(NAME [soildynamics]2d_semi-circle_transient_traction COMMAND Soildynamics
    inputs/transient-traction.arc)
  add_test(NAME [soildynamics]2d_bar_elastodynamics COMMAND Soildynamics
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/Soildynamics.arc)
  add_test(NAME [soildynamics]2d_bar_dirichlet_rowelm COMMAND Soildynamics
    -A,//fem/boundary-conditions/dirichlet[1]/enforce-Dirichlet-method=RowElimination
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/Soildynamics.arc)
  add_test(NAME [soildynamics]parse_and_exit COMMAND Soildynamics
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    -A,//fem/cross-validation=false
    inputs/Soildynamics.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [soildynamics]2d_semi-circle_const_traction_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics inputs/constant-traction.arc)
    add_test(NAME [soildynamics]2d_semi-circle_transient_traction_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics inputs/transient-traction.arc)
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
  -A,//fem/linear-system/@name=HypreLinearSystem
  -A,//fem/linear-system/solver=cg
  -A,//fem/linear-system/atol=1e-9
  -A,//fem/linear-system/rtol=0)

  add_test(NAME [soildynamics]2d_double-couple_paraxial COMMAND Soildynamics
    inputs/double-couple.paraxial.arc)

  add_test(NAME [soildynamics]2d_double-couple_paraxial_body-force COMMAND Soildynamics
    -A,//fem/linear-system/solver-backend=hypre
    inputs/double-couple.paraxial.body-force.arc)

  add_test(NAME [soildynamics]2d_double-couple_paraxial_soil COMMAND Soildynamics
    inputs/double-couple.paraxial.soil.arc)

  add_test(NAME [soildynamics]3d_double-couple_paraxial_soil COMMAND Soildynamics
    -A,//fem/linear-system/solver-backend=hypre
    inputs/3d.double-couple.paraxial.soil.arc)

  add_test(NAME [soildynamics]3d_all-conditions COMMAND Soildynamics
    inputs/3d.all-conditions.arc)

  arcanefem_add_gpu_test(NAME [soildynamics]2d_double-couple_paraxial_body-force_bsr COMMAND  ./Soildynamics ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/double-couple.paraxial.body-force.arc)

  arcanefem_add_gpu_test(NAME [soildynamics]2d_double-couple_paraxial_body-force_af-bsr COMMAND  ./Soildynamics ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/double-couple.paraxial.body-force.arc)

  arcanefem_add_gpu_test(NAME [soildynamics]3d_double-couple_paraxial_soil_bsr COMMAND  ./Soildynamics ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/3d.double-couple.paraxial.soil.arc)

  arcanefem_add_gpu_test(NAME [soildynamics]3d_double-couple_paraxial_soil_af-bsr COMMAND  ./Soildynamics ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/3d.double-couple.paraxial.soil.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [soildynamics]2d_double-couple_paraxial_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics inputs/double-couple.paraxial.arc)
    add_test(NAME [soildynamics]2d_double-couple_paraxial_soil_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics inputs/double-couple.paraxial.soil.arc)
    add_test(NAME [soildynamics]3d_double-couple_paraxial_soil_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics
      -A,//fem/linear-system/solver-backend=hypre
      inputs/3d.double-couple.paraxial.soil.arc)
    add_test(NAME [soildynamics]3d_all-conditions_2p  COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics inputs/3d.all-conditions.arc)
  endif()
endif()
