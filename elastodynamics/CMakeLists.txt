add_executable(Elastodynamics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elastodynamics)
target_include_directories(Elastodynamics PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Elastodynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elastodynamics.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elastodynamics.pointBC.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elastodynamics.damping.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elastodynamics.dirichlet.traction.bodyforce.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elastodynamics.Galpha.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Elastodynamics.transient-traction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(traction_bar_test_1.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/semi-circle.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elastodynamics PUBLIC FemUtils)

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [elastodynamics] COMMAND Elastodynamics Test.Elastodynamics.arc)
  add_test(NAME [elastodynamics]transient_traction COMMAND Elastodynamics Test.Elastodynamics.transient-traction.arc)
  add_test(NAME [elastodynamics]Dirichlet_pointBc COMMAND Elastodynamics Test.Elastodynamics.pointBC.arc)
  add_test(NAME [elastodynamics]constant_traction_and_damping COMMAND Elastodynamics Test.Elastodynamics.damping.arc)
  add_test(NAME [elastodynamics]time-discretization_Galpha COMMAND Elastodynamics Test.Elastodynamics.Galpha.arc)
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [elastodynamics]Dirichlet_traction_bodyforce COMMAND Elastodynamics Test.Elastodynamics.dirichlet.traction.bodyforce.arc)
endif()
