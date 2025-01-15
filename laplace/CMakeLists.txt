set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Laplace
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Laplace)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Laplace)
target_include_directories(Laplace PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Laplace.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.ring.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.L-shape.3D.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.PointDirichlet.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.PointDirichlet-refined.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/ring.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/L-shape-3D.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Laplace PUBLIC FemUtils)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plancher.geo -format msh41 -o plancher.10k.msh -bin
    DEPENDS ${MSH_DIR}/plancher.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files_laplace DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh)
  add_dependencies(Laplace gmsh_files_laplace)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [laplace]2D_Dirichlet COMMAND Laplace Test.ring.arc)
endif()

add_test(NAME [laplace]3D_Dirichlet COMMAND Laplace Test.L-shape.3D.arc)
add_test(NAME [laplace]pointDirichlet COMMAND Laplace Test.PointDirichlet.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [laplace]2D_Dirichlet_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace Test.ring.arc)
  add_test(NAME [laplace]3D_Dirichlet_3D_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace Test.L-shape.3D.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME [laplace]pointDirichlet_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace Test.PointDirichlet-refined.arc)
  endif()
endif()
