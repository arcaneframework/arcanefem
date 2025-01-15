add_executable(Elastodynamics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elastodynamics)
target_include_directories(Elastodynamics PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Elastodynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elastodynamics PUBLIC FemUtils)

# Copy the data files in the binary directory
file(COPY "data" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  bar_dynamic.msh
  semi-circle.msh 
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()


enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [elastodynamics] COMMAND Elastodynamics inputs/bar.arc)
  add_test(NAME [elastodynamics]transient_traction COMMAND Elastodynamics inputs/bar.transient-traction.arc)
  add_test(NAME [elastodynamics]Dirichlet_pointBc COMMAND Elastodynamics inputs/semi-circle.pointBC.arc)
  add_test(NAME [elastodynamics]constant_traction_and_damping COMMAND Elastodynamics inputs/bar.damping.arc)
  add_test(NAME [elastodynamics]time-discretization_Galpha COMMAND Elastodynamics inputs/bar.Galpha.arc)
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [elastodynamics]Dirichlet_traction_bodyforce COMMAND Elastodynamics inputs/bar.dirichlet.traction.bodyforce.arc)
endif()
