add_executable(Soildynamics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Soildynamics)
target_include_directories(Soildynamics PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Soildynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Soildynamics.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.transient-traction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.constant-traction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.double-couple.paraxial.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.double-couple.paraxial.soil.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/semi-circle-soil.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/square_double-couple.geo ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/soil_2d.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Soildynamics PUBLIC FemUtils)

file(COPY "data/"  DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [soildynamics]soildynamics COMMAND Soildynamics Test.Soildynamics.arc)
  add_test(NAME [soildynamics]soildynamics_const_traction COMMAND Soildynamics Test.constant-traction.arc)
  add_test(NAME [soildynamics]soildynamics_transient_traction COMMAND Soildynamics Test.transient-traction.arc)
  add_test(NAME [soildynamics]soildynamics_dc_paraxial COMMAND Soildynamics Test.double-couple.paraxial.soil.arc)
endif()