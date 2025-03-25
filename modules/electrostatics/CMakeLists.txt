add_executable(Electrostatics
  FemModule.cc
  main.cc
  Fem_axl.h
)

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
  box-rod-circle.msh
  interdigital_capacitor.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [electrostatics]box-rods COMMAND Electrostatics inputs/box-rods.arc)
  add_test(NAME [electrostatics]capacitor COMMAND Electrostatics inputs/Capacitor.arc)
  add_test(NAME [electrostatics]rod-circle COMMAND Electrostatics inputs/rod-circle.arc)
  add_test(NAME [electrostatics]box-rods_petsc_monitor COMMAND Electrostatics
    -A,//fem/petsc-flags=-ksp_monitor
    inputs/box-rods.arc)
  add_test(NAME [electrostatics]2D_parse_and_exit COMMAND Electrostatics
    -A,//fem/assemble-linear-system=false
    -A,//fem/solve-linear-system=false
    -A,//fem/cross-validation=false
    inputs/Capacitor.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [electrostatics]box-rods_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Electrostatics inputs/box-rods.arc)
  endif()
endif()
