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

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  bar_dynamic.msh
  bar_dynamic_3D.msh
  semi-circle.msh 
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()


enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [elastodynamics]2D_bar COMMAND Elastodynamics inputs/bar.arc)
  add_test(NAME [elastodynamics]2D_bar_RowElimination COMMAND Elastodynamics
    -A,//fem/enforce-Dirichlet-method=RowElimination
    inputs/bar.arc)
  add_test(NAME [elastodynamics]2D_bar_RowColumnElimination COMMAND Elastodynamics
    -A,//fem/enforce-Dirichlet-method=RowColumnElimination
    inputs/bar.arc)
  add_test(NAME [elastodynamics]3D_bar COMMAND Elastodynamics inputs/bar.3D.arc)
  add_test(NAME [elastodynamics]2D_transient_traction COMMAND Elastodynamics inputs/bar.transient-traction.arc)
  add_test(NAME [elastodynamics]3D_transient_traction COMMAND Elastodynamics inputs/bar.3D.transient-traction.arc)
  add_test(NAME [elastodynamics]2D_constant_traction_and_damping COMMAND Elastodynamics inputs/bar.damping.arc)
  add_test(NAME [elastodynamics]2D_time-discretization_Galpha COMMAND Elastodynamics inputs/bar.Galpha.arc)
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [elastodynamics]Dirichlet_traction_bodyforce COMMAND Elastodynamics inputs/bar.dirichlet.traction.bodyforce.arc)
  add_test(NAME [elastodynamics]Dirichlet_pointBc_Penalty COMMAND Elastodynamics
    -A,//fem/enforce-Dirichlet-method=Penalty
    inputs/semi-circle.pointBC.arc)
  add_test(NAME [elastodynamics]2D_pointBc_RowElimination COMMAND Elastodynamics
    -A,//fem/enforce-Dirichlet-method=RowElimination
    inputs/semi-circle.pointBC.arc)
  add_test(NAME [elastodynamics]2D_pointBc_RowColumnElimination COMMAND Elastodynamics
    -A,//fem/enforce-Dirichlet-method=RowColumnElimination
    inputs/semi-circle.pointBC.arc)
endif()
