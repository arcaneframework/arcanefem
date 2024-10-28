add_executable(Passmo
  TypesElastodynamic.h
  Elastodynamic_axl.h
  ElastodynamicModule.h
  ElastodynamicModule.cc
  main.cc
  analytic_func.cc
  analytical_func.h)

arcane_generate_axl(Elastodynamic)
arcane_add_arcane_libraries_to_target(Passmo)

target_include_directories(Passmo PUBLIC . ../fem1 ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Passmo.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_hybrid_quad_tria.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Passmo PUBLIC FemUtils)

file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "inputs/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()
if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [passmo]bar2D_T3 COMMAND Passmo bar2d_tria.arc)
  add_test(NAME [passmo]bar2D_Q4 COMMAND Passmo bar2d_quad.arc)
  add_test(NAME [passmo]bar2D_T3-Q4 COMMAND Passmo bar2d_tria_quad.arc)
endif ()

if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [passmo]bar2D_T3_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmo bar2d_tria.arc)
  add_test(NAME [passmo]bar2D_T3-Q4_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmo bar2d_tria_quad.arc)
endif()
