set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Poisson_dg
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Poisson_dg)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Poisson_dg)
target_include_directories(Poisson_dg PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Poisson_dg.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Poisson_dg PUBLIC FemUtils)

# Copy the check files in the binary directory
# file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
    circle_cut.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

add_test(NAME [Poisson_dg]circle_homogenous_tria COMMAND Poisson_dg
  inputs/circle.2D.dg.arc)