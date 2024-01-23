add_executable(Elasticity
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elasticity)
target_include_directories(Elasticity PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Elasticity.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elasticity.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elasticity.traction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elasticity.PointDirichlet.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elasticity.DirichletViaRowElimination.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elasticity.DirichletViaRowColumnElimination.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elasticity PUBLIC FemUtils)

enable_testing()
add_test(NAME [elasticity] COMMAND Elasticity Test.Elasticity.arc)
add_test(NAME [elasticity]Dirichlet_traction COMMAND Elasticity Test.Elasticity.traction.arc)
add_test(NAME [elasticity]Dirichlet_pointBC COMMAND Elasticity Test.Elasticity.PointDirichlet.arc)
add_test(NAME [elasticity]Dirichlet_via_RowElimination COMMAND Elasticity Test.Elasticity.DirichletViaRowElimination.arc)
add_test(NAME [elasticity]Dirichlet_via_RowColElimination COMMAND Elasticity Test.Elasticity.DirichletViaRowColumnElimination.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [elasticity]parallel_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity Test.Elasticity.arc)
  add_test(NAME [elasticity]parallel_Dirichlet_RowElimination_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity Test.Elasticity.DirichletViaRowElimination.arc)
  add_test(NAME [elasticity]parallel_Dirichlet_RowColElimination_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Elasticity Test.Elasticity.DirichletViaRowColumnElimination.arc)
endif()
