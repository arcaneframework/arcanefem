add_executable(heat
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(heat)
target_include_directories(heat PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Heat.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  plate.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

target_link_libraries(heat PUBLIC FemUtils)

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/meshes/plate.fine.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plate.fine.geo -format msh41 -o meshes/plate.fine.msh -bin
    DEPENDS ${MSH_DIR}/plate.fine.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files_heat DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/meshes/plate.fine.msh)
  add_dependencies(heat gmsh_files_heat)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()


if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [heat]conduction COMMAND heat inputs/conduction.arc)
  add_test(NAME [heat]conduction_RowElimination_Dirichlet COMMAND heat inputs/conduction.DirichletViaRowElimination.arc)
  add_test(NAME [heat]conduction_RowColElimination_Dirichlet COMMAND heat inputs/conduction.DirichletViaRowColumnElimination.arc)
  add_test(NAME [heat]conduction_convection COMMAND heat inputs/conduction.convection.arc)
endif()


# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [heat]conduction_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.arc)
  add_test(NAME [heat]conduction_RowElimination_Dirichlet_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.DirichletViaRowElimination.arc)
  add_test(NAME [heat]conduction_RowColElimination_Dirichlet_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.DirichletViaRowColumnElimination.arc)
  add_test(NAME [heat]conduction_convection_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.convection.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME [heat]conduction_fine_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.fine.arc)
    add_test(NAME [heat]conduction_convection_fine_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./heat inputs/conduction.convection.fine.arc)
  endif()
endif()
