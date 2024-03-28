add_executable(Passmo
  TypesElastodynamic.h
  Elastodynamic_axl.h
  ElastodynamicModule.h
  ElastodynamicModule.cc
  main.cc
  utilFEM.h
  utilFEM.cc
  Integer3std.h
  analytic_func.cc
  analytical_func.h
  GaussDoFsOnCells.h
  GaussDoFsOnCells.cc
  GaussQuadrature.h
  GaussQuadrature.cc )

arcane_generate_axl(Elastodynamic)
arcane_add_arcane_libraries_to_target(Passmo)

target_include_directories(Passmo PUBLIC . ../fem1 ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Passmo.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(Test.Passmo.transient-traction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(Test.Passmo.constant-traction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(semi-circle-soil-traction.txt ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
#configure_file(${MSH_DIR}/semi-circle-soil.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Passmo PUBLIC FemUtils)


#enable_testing()
#add_test(NAME [passmo]passmo_const_traction COMMAND Passmo Test.Passmo.constant-traction.arc)
#add_test(NAME [passmo]passmo_transient_traction COMMAND Passmo Test.Passmo.transient-traction.arc)
