add_executable(Soildynamics
  FemModule.cc
  main.cc
  Fem_axl.h
)

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
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [soildynamics]soildynamics COMMAND Soildynamics inputs/Soildynamics.arc)
  add_test(NAME [soildynamics]soildynamics_const_traction COMMAND Soildynamics inputs/constant-traction.arc)
  add_test(NAME [soildynamics]soildynamics_const_traction_pointbc COMMAND Soildynamics inputs/constant-traction.pointbc.arc)
  add_test(NAME [soildynamics]soildynamics_transient_traction COMMAND Soildynamics inputs/transient-traction.arc)
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [soildynamics]soildynamics_dc_paraxial COMMAND Soildynamics inputs/double-couple.paraxial.arc)
  add_test(NAME [soildynamics]soildynamics_dc_paraxial_body-force COMMAND Soildynamics inputs/double-couple.paraxial.body-force.arc)
  add_test(NAME [soildynamics]soildynamics_dc_paraxial_soil COMMAND Soildynamics inputs/double-couple.paraxial.soil.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [soildynamics]soildynamics_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics inputs/Soildynamics.arc)
    add_test(NAME [soildynamics]soildynamics_dc_paraxial_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Soildynamics inputs/double-couple.paraxial.arc)
  endif()
endif()
