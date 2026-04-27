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
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MSH_FILES
    circle_cut.msh
)
set(MED_FILES
    circle_cut-tria.med
    circle_cut-quad-tria.med
    circle_cut-quad.med
)
foreach(MSH_FILES IN LISTS MSH_FILES)
    file(COPY ${MSH_DIR}/${MSH_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()
foreach(MED_FILES IN LISTS MED_FILES)
    file(COPY ${MED_DIR}/${MED_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()


enable_testing()

add_test(NAME [Poisson_dg]circle_homogenous_tria COMMAND Poisson_dg
  -A,//fem/solution-comparison-file=check/circle_2D.txt
  inputs/circle.2D.arc)

if(MEDFile_FOUND)
  message(STATUS "MEDFile found: ${MEDFile_FOUND}")
  add_test(NAME [Poisson_dg]circle_homogenous_tria_med COMMAND Poisson_dg
    -A,//meshes/mesh/filename=meshes/circle_cut-tria.med
#    -A,//fem/solution-comparison-file=check/circle_2D.txt
    inputs/circle.2D.arc)
  add_test(NAME [Poisson_dg]circle_homogenous_quad_tria_med COMMAND Poisson_dg
    -A,//meshes/mesh/filename=meshes/circle_cut-quad-tria.med
    inputs/circle.2D.arc)
  add_test(NAME [Poisson_dg]circle_homogenous_quad_med COMMAND Poisson_dg
    -A,//meshes/mesh/filename=meshes/circle_cut-quad.med
    inputs/circle.2D.arc)
endif()