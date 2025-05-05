add_executable(Passmonl
  TypesNLdynamic.h
  NLDynamic_axl.h
  NLDynamicModule.h
  NLDynamicModule.cc
  main.cc
  analytic_func.cc
  analytical_func.h
        utilFEM.h)

arcane_generate_axl(NLDynamic)
arcane_add_arcane_libraries_to_target(Passmonl)

target_include_directories(Passmonl PUBLIC . ../fem1 ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(Passmonl PUBLIC FemUtils)
configure_file(Passmonl.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

configure_file(${MSH_DIR}/bar_hybrid_quad_tria.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic_quad.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic_3D.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic_3Dhexa.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "inputs/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#enable_testing()
#if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
#  add_test(NAME [passmonl]bar2D_T3 COMMAND Passmonl bar2d_tria.arc)
#  add_test(NAME [passmonl]bar2D_Q4 COMMAND Passmonl bar2d_quad.arc)
#  add_test(NAME [passmonl]bar2D_T3-Q4 COMMAND Passmonl bar2d_tria_quad.arc)
#  add_test(NAME [passmonl]bar3D_T4 COMMAND Passmonl bar3d_tetra.arc)
#  add_test(NAME [passmonl]bar3D_H8 COMMAND Passmonl bar3d_hexa.arc)
#endif ()

#if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
#  add_test(NAME [passmonl]bar2D_T3 COMMAND Passmonl bar2d_tria.arc)
#  add_test(NAME [passmonl]bar2D_Q4 COMMAND Passmonl bar2d_quad.arc)
#  add_test(NAME [passmonl]bar2D_T3-Q4 COMMAND Passmonl bar2d_tria_quad.arc)
#  add_test(NAME [passmonl]bar3D_T4 COMMAND Passmonl bar3d_tetra.arc)
#  add_test(NAME [passmonl]bar3D_H8 COMMAND Passmonl bar3d_hexa.arc)
#endif ()

if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [passmonl]bar2D_T3_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmonl bar2d_tria.arc)
  add_test(NAME [passmonl]bar2D_Q4_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmonl bar2d_quad.arc)
  add_test(NAME [passmonl]bar2D_T3-Q4_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmonl bar2d_tria_quad.arc)
  add_test(NAME [passmonl]bar3D_T4_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmonl bar3d_tetra.arc)
  add_test(NAME [passmonl]bar3D_H8_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmonl bar3d_hexa.arc)
endif()
