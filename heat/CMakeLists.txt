add_executable(heat
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(heat)
target_include_directories(heat PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Heat.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.DirichletViaRowElimination.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.DirichletViaRowColumnElimination.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.convection.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.fine.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.convection.fine.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/plate.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(heat PUBLIC FemUtils)

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/plate.fine.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plate.fine.geo -format msh41 -o plate.fine.msh
    DEPENDS ${MSH_DIR}/plate.fine.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files_heat DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/plate.fine.msh)
  add_dependencies(heat gmsh_files_heat)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()


if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [heat]conduction COMMAND heat Test.conduction.arc)
  add_test(NAME [heat]conduction_RowElimination_Dirichlet COMMAND heat Test.conduction.DirichletViaRowElimination.arc)
  add_test(NAME [heat]conduction_RowColElimination_Dirichlet COMMAND heat Test.conduction.DirichletViaRowColumnElimination.arc)
  add_test(NAME [heat]conduction_convection COMMAND heat Test.conduction.convection.arc)
endif()


# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [heat]conduction_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.arc)
  add_test(NAME [heat]conduction_RowElimination_Dirichlet_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.DirichletViaRowElimination.arc)
  add_test(NAME [heat]conduction_RowColElimination_Dirichlet_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.DirichletViaRowColumnElimination.arc)
  add_test(NAME [heat]conduction_convection_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.convection.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME [heat]conduction_fine_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.fine.arc)
    add_test(NAME [heat]conduction_convection_fine_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat Test.conduction.convection.fine.arc)
  endif()
endif()
