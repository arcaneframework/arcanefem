set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Electrostatics
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Electrostatics)
target_include_directories(Electrostatics PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Electrostatics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Electrostatics PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  box-rods.msh
  box-rods.quad.msh
  box-rod-circle.msh
  interdigital_capacitor.msh
  truncated_cube.msh
  truncated_cube.hexa.msh
  spheres_in_sphere.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [electrostatics]box-rods COMMAND Electrostatics inputs/box-rods.arc)
  add_test(NAME [electrostatics]box-rods_quad COMMAND Electrostatics inputs/box-rods.quad.arc)
  add_test(NAME [electrostatics]capacitor COMMAND Electrostatics inputs/Capacitor.arc)
  add_test(NAME [electrostatics]rod-circle COMMAND Electrostatics inputs/rod-circle.arc)
  add_test(NAME [electrostatics]spheres COMMAND Electrostatics inputs/spheres.arc)
  add_test(NAME [electrostatics]truncated-cube_hexa COMMAND Electrostatics inputs/truncated_cube.hexa.arc)
  add_test(NAME [electrostatics]box-rods_petsc_monitor COMMAND Electrostatics
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/box-rods.arc)
  add_test(NAME [electrostatics]2D_parse_and_exit COMMAND Electrostatics
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    -A,//fem/cross-validation=false
    inputs/Capacitor.arc)

  arcanefem_add_gpu_test(NAME [electrostatics]box-rods_bsr COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=BSR
    inputs/box-rods.arc)
  arcanefem_add_gpu_test(NAME [electrostatics]box-rods_af-bsr COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/box-rods.arc)
  arcanefem_add_gpu_test(NAME [electrostatics]truncated_cube_bsr COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=BSR
    inputs/truncated_cube.arc)
  arcanefem_add_gpu_test(NAME [electrostatics]spheres_bsr COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=BSR
    inputs/spheres.arc)
  arcanefem_add_gpu_test(NAME [electrostatics]spheres_af-bsr COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/spheres.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [electrostatics]box-rods_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Electrostatics inputs/box-rods.arc)
    add_test(NAME [electrostatics]box-rods_quad_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Electrostatics inputs/box-rods.quad.arc)
    add_test(NAME [electrostatics]truncated-cube_hexa_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Electrostatics inputs/truncated_cube.hexa.arc)
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)
  arcanefem_add_gpu_test(NAME [electrostatics]rod-circle_bsr_hypreDirect COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/rod-circle.arc)
  arcanefem_add_gpu_test(NAME [electrostatics]rod-circle_af-bsr_hypreDirect COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/rod-circle.arc)
  arcanefem_add_gpu_test(NAME [electrostatics]spheres_af-bsr_hypreDirect COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/spheres.arc)
  arcanefem_add_gpu_test(NAME [electrostatics]truncated_cube_bsr_hypreDirect COMMAND ./Electrostatics ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/truncated_cube.arc)
endif()
