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
configure_file(${MSH_DIR}/bar_hybrid_quad_tria.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Passmo.bar2d_hybridpar.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.bar2d_tria.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Passmo PUBLIC FemUtils)


enable_testing()
add_test(NAME [passmo]2D_hybrid_T3-Q4 COMMAND Passmo Test.Passmo.bar2d_hybridpar.arc)
add_test(NAME [passmo]2D_T3 COMMAND Passmo Test.bar2d_tria.arc)

if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [passmo]2D_T3_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmo Test.bar2d_tria.arc)
endif()
