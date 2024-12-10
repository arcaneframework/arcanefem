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
configure_file(Test.bar.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.bar.traction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.bar.PointDirichlet.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.bar.DirichletViaRowElimination.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.bar.DirichletViaRowColumnElimination.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elasticity PUBLIC FemUtils)

enable_testing()

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [elasticity] COMMAND Elasticity Test.bar.arc)
  add_test(NAME [elasticity]Dirichlet_traction COMMAND Elasticity Test.bar.traction.arc)
  add_test(NAME [elasticity]Dirichlet_pointBC COMMAND Elasticity Test.bar.PointDirichlet.arc)
  add_test(NAME [elasticity]Dirichlet_via_RowElimination COMMAND Elasticity Test.bar.DirichletViaRowElimination.arc)
  add_test(NAME [elasticity]Dirichlet_via_RowColElimination COMMAND Elasticity Test.bar.DirichletViaRowColumnElimination.arc)
endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  # Temporarely remove this test because there is a difference on node 37
  add_test(NAME [elasticity]parallel_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity Test.bar.arc)
  add_test(NAME [elasticity]parallel_Dirichlet_RowElimination_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity Test.bar.DirichletViaRowElimination.arc)
  add_test(NAME [elasticity]parallel_Dirichlet_RowColElimination_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity Test.bar.DirichletViaRowColumnElimination.arc)
endif()
