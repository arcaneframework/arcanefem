add_executable(Poisson
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Poisson)
target_include_directories(Poisson PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Poisson.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Poisson PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
    sphere_cut.msh
    circle_cut.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()


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
