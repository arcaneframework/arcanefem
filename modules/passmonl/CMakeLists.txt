add_executable(Passmonl
  TypesNLdynamic.h
  NLDynamic_axl.h
  NLDynamicModule.h
  NLDynamicModule.cc
  main.cc
  utilFEM.h
  utilFEM.cc
  LawDispatcher.h
  LawDispatcher.cc
  analytic_func.cc
  analytical_func.h
  hooke.cc
  druckerp.cc
#        BoundaryConditions.cc
#        GaussItem.h
)

arcane_generate_axl(NLDynamic)
arcane_add_arcane_libraries_to_target(Passmonl)

target_include_directories(Passmonl PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(Passmonl PUBLIC FemUtils)
configure_file(Passmonl.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

configure_file(${MSH_DIR}/bar.two_trias.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar4Q.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/bar_dynamic.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file(${MSH_DIR}/bar_dynamic_quad.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file(${MSH_DIR}/bar_dynamic_3D.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# configure_file(${MSH_DIR}/bar_dynamic_3Dhexa.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# Copy the tests and the inputs folder containing the arc files
#file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY "inputs/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()
if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [passmonl]bar2D_T3_petsc COMMAND Passmonl bar2d_2T_lin_2steps.arc)
  add_test(NAME [passmonl]bar2D_Q4_dp_petsc COMMAND Passmonl bar2d_4Q_dp.arc)
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [passmonl]bar2D_T3_hypre COMMAND Passmonl bar2d_2T_lin_2steps.arc)
  add_test(NAME [passmonl]bar2D_Q4_dp_hypre COMMAND Passmonl bar2d_4Q_dp.arc)
endif ()

if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [passmonl]bar2D_T3_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmonl
    -A,//n-l-dynamic/linear-system/@name=AlephLinearSystem
    bar2d_2T_lin_2steps.arc)
  # add_test(NAME [passmonl]bar2D_Q4_dp_2pu COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Passmonl
  #   -A,//n-l-dynamic/linear-system/@name=AlephLinearSystem
  #   bar2d_4Q_dp.arc)
endif()
#enable_testing()
#set(MODULE_NAME "passmonl")

#if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
#  arcanefem_add_seq_test(NAME bar2D_T3 COMMAND ./Passmonl ARGS bar2d_tria.arc)
#  arcanefem_add_seq_test(NAME bar2D_Q4 COMMAND ./Passmonl ARGS bar2d_quad.arc)
#  arcanefem_add_seq_test(NAME bar2D_T3-Q4 COMMAND ./Passmonl ARGS bar2d_tria_quad.arc)
#  arcanefem_add_seq_test(NAME bar3D_T4 COMMAND ./Passmonl ARGS bar3d_tetra.arc)
#  arcanefem_add_seq_test(NAME bar3D_H8 COMMAND ./Passmonl ARGS bar3d_hexa.arc)
#endif ()

#if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
#  arcanefem_add_seq_test(NAME bar2D_T3 COMMAND ./Passmonl ARGS bar2d_tria.arc)
#  arcanefem_add_seq_test(NAME bar2D_Q4 COMMAND ./Passmonl ARGS bar2d_quad.arc)
#  arcanefem_add_seq_test(NAME bar2D_T3-Q4 COMMAND ./Passmonl ARGS bar2d_tria_quad.arc)
#  arcanefem_add_seq_test(NAME bar3D_T4 COMMAND ./Passmonl ARGS bar3d_tetra.arc)
#  arcanefem_add_seq_test(NAME bar3D_H8 COMMAND ./Passmonl ARGS bar3d_hexa.arc)
#endif ()

#if(FEMUTILS_HAS_PARALLEL_SOLVER)
#   arcanefem_add_mpi_test(NAME bar2D_T3 COMMAND ./Passmonl ARGS bar2d_tria.arc)
#   arcanefem_add_mpi_test(NAME bar2D_Q4 COMMAND ./Passmonl ARGS bar2d_quad.arc)
#   arcanefem_add_mpi_test(NAME bar2D_T3-Q4 COMMAND ./Passmonl ARGS bar2d_tria_quad.arc)
#   arcanefem_add_mpi_test(NAME bar3D_T4 COMMAND ./Passmonl ARGS bar3d_tetra.arc)
#   arcanefem_add_mpi_test(NAME bar3D_H8 COMMAND ./Passmonl ARGS bar3d_hexa.arc)
# endif()
