add_executable(Poisson
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Poisson)
target_include_directories(Poisson PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Poisson.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)


# Copy the mesh files
configure_file(${MSH_DIR}/sphere_cut.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/circle_cut.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Poisson PUBLIC FemUtils)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [poisson]2D COMMAND Poisson inputs/circle.2D.arc)
  add_test(NAME [poisson]3D COMMAND Poisson inputs/sphere.3D.arc)
endif()



# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [poisson]2D_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson inputs/circle.2D.arc)
  add_test(NAME [poisson]3D_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson inputs/sphere.3D.arc)
endif()
