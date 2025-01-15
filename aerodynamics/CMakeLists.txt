add_executable(aerodynamics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(aerodynamics)
target_include_directories(aerodynamics PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(aerodynamics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(aerodynamics PUBLIC FemUtils)

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  NACA0012.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [aerodynamics]Joukowski COMMAND aerodynamics inputs/Joukowski.arc)
endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [aerodynamics]Joukowski_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./aerodynamics inputs/Joukowski.arc)
endif()
