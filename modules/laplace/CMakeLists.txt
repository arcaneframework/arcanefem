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

target_link_libraries(Laplace PUBLIC FemUtils)

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  ring.msh
  plancher.msh
  L-shape-3D.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/meshes/plancher.10k.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plancher.geo -format msh41 -o meshes/plancher.10k.msh -bin
    DEPENDS ${MSH_DIR}/plancher.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(gmsh_files_laplace DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/meshes/plancher.10k.msh)
  add_dependencies(Laplace gmsh_files_laplace)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [laplace]2D_Dirichlet COMMAND Laplace inputs/ring.arc)
  add_test(NAME [laplace]3D_Dirichlet COMMAND Laplace inputs/L-shape.3D.arc)
endif()

add_test(NAME [laplace]pointDirichlet COMMAND Laplace inputs/PointDirichlet.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [laplace]pointDirichlet_bsr_hypreDirect COMMAND Laplace inputs/PointDirichlet.bsr.hypreDirect.arc)
  arcanefem_add_gpu_test(NAME [laplace]pointDirichlet_bsr_hypreDirect_gpu COMMAND ./Laplace ARGS inputs/PointDirichlet.bsr.atomicFree.hypreDirect.arc)
  add_test(NAME [laplace]pointDirichlet_bsr_atomicFree_hypreDirect COMMAND Laplace inputs/PointDirichlet.bsr.hypreDirect.arc)
  arcanefem_add_gpu_test(NAME [laplace]pointDirichlet_bsr_atomicFree_hypreDirect_gpu COMMAND ./Laplace ARGS inputs/PointDirichlet.bsr.atomicFree.hypreDirect.arc)
  endif()

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [laplace]2D_Dirichlet_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace inputs/ring.arc)
  add_test(NAME [laplace]3D_Dirichlet_3D_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace inputs/L-shape.3D.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME [laplace]pointDirichlet_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Laplace inputs/PointDirichlet-refined.arc)
  endif()
endif()
