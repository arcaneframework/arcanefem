add_executable(Acoustics
  FemModule.h
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Acoustics)

target_include_directories(Acoustics PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Acoustics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  sub.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

target_link_libraries(Acoustics PUBLIC FemUtils)

enable_testing()

add_test(NAME [Acoustics]2D_submarine COMMAND Acoustics inputs/sub.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [Acoustics]2D_submarine_hypre COMMAND Acoustics inputs/sub.hypre.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [Acoustics]2D_submarine_hypre_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Acoustics inputs/sub.hypre.arc)
  endif()
endif()
