set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Elasticity
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Elasticity)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elasticity)

target_include_directories(Elasticity PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})

configure_file(Elasticity.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elasticity PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
    bar.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()


enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [elasticity]Dirichlet COMMAND Elasticity inputs/bar.2D.arc)
  add_test(NAME [elasticity]Dirichlet_traction COMMAND Elasticity inputs/bar.2D.traction.arc)
  add_test(NAME [elasticity]Dirichlet_pointBC COMMAND Elasticity inputs/bar.2D.PointDirichlet.arc)
  add_test(NAME [elasticity]Dirichlet_pointBC_bsr COMMAND Elasticity inputs/bar.2D.PointDirichlet.bsr.arc)
  add_test(NAME [elasticity]Dirichlet_pointBC_via_RowElimination COMMAND Elasticity inputs/bar.2D.PointDirichlet.DirichletViaRowElimination.arc)
  add_test(NAME [elasticity]Dirichlet_pointBC_via_RowColElimination COMMAND Elasticity inputs/bar.2D.PointDirichlet.DirichletViaRowColumnElimination.arc)
  add_test(NAME [elasticity]Dirichlet_via_RowElimination COMMAND Elasticity inputs/bar.2D.DirichletViaRowElimination.arc)
  add_test(NAME [elasticity]Dirichlet_via_RowColElimination COMMAND Elasticity inputs/bar.2D.DirichletViaRowColumnElimination.arc)
  add_test(NAME [elasticity]Dirichlet_traction_bsr COMMAND Elasticity inputs/bar.2D.traction.bsr.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_bsr_gpu COMMAND ./Elasticity ARGS inputs/bar.2D.traction.bsr.arc)
  add_test(NAME [elasticity]Dirichlet_traction_Bodyforce COMMAND Elasticity inputs/bar.2D.traction.bodyforce.arc)
  add_test(NAME [elasticity]Dirichlet_traction_Bodyforce_bsr COMMAND Elasticity inputs/bar.2D.traction.bodyforce.bsr.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_Bodyforce_bsr_gpu COMMAND ./Elasticity ARGS inputs/bar.2D.traction.bodyforce.bsr.arc)
  add_test(NAME [elasticity]Dirichlet_traction_Bodyforce_bsr_atomic_free COMMAND Elasticity inputs/bar.2D.traction.bodyforce.bsr.atomic-free.arc)
  arcanefem_add_gpu_test(NAME [elasticity]Dirichlet_traction_Bodyforce_bsr_gpu_atomic_free COMMAND ./Elasticity ARGS inputs/bar.2D.traction.bodyforce.bsr.atomic-free.arc)
endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  # Temporarely remove this test because there is a difference on node 37
  add_test(NAME [elasticity]Dirichlet_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity inputs/bar.2D.arc)
  add_test(NAME [elasticity]Dirichlet_via_RowElimination_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity inputs/bar.2D.DirichletViaRowElimination.arc)
  add_test(NAME [elasticity]Dirichlet_via_RowColElimination_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity inputs/bar.2D.DirichletViaRowColumnElimination.arc)
endif()
