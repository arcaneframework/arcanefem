set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Poisson
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Poisson)

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
    plancher.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()


enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [poisson]2D COMMAND Poisson inputs/circle.2D.arc)
  add_test(NAME [poisson]2D_neumann COMMAND Poisson inputs/circle.neumann.2D.arc)

  add_test(NAME [poisson]2D_bsr COMMAND Poisson inputs/circle.2D.bsr.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_bsr_gpu COMMAND Poisson ARGS inputs/circle.2D.bsr.arc)

  add_test(NAME [poisson]2D_bsr_atomicFree COMMAND Poisson inputs/circle.2D.bsr.atomicFree.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_bsr_atomicFree_gpu COMMAND Poisson ARGS inputs/circle.2D.bsr.atomicFree.arc)

  add_test(NAME [poisson]3D COMMAND Poisson inputs/sphere.3D.arc)
  add_test(NAME [poisson]3D_neumann COMMAND Poisson inputs/sphere.neumann.3D.arc)

  add_test(NAME [poisson]3D_bsr COMMAND Poisson inputs/sphere.3D.bsr.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_bsr_gpu COMMAND Poisson ARGS inputs/sphere.3D.bsr.arc)

  add_test(NAME [poisson]3D_bsr_atomicFree COMMAND Poisson inputs/sphere.3D.bsr.atomicFree.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_bsr_atomicFree_gpu COMMAND Poisson ARGS inputs/sphere.3D.bsr.atomicFree.arc)
endif()

if (FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [poisson]2D_pntDirichlet_bsr_hypre COMMAND Poisson inputs/perforatedSquare.pointDirichlet.2D.bsr.hypre.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_pntDirichlet_bsr_hypre_gpu COMMAND Poisson ARGS inputs/perforatedSquare.pointDirichlet.2D.bsr.hypre.arc)

  add_test(NAME [poisson]2D_bsr_hypre COMMAND Poisson inputs/circle.2D.bsr.hypre.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_bsr_hypre_gpu COMMAND Poisson ARGS inputs/circle.2D.bsr.hypre.arc)

  add_test(NAME [poisson]2D_neumann_bsr_hypre COMMAND Poisson inputs/circle.neumann.2D.bsr.hypre.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_neumann_bsr_hypre_gpu COMMAND Poisson ARGS inputs/circle.neumann.2D.bsr.hypre.arc)

  add_test(NAME [poisson]2D_bsr_atomicFree_hypre COMMAND Poisson inputs/circle.2D.bsr.atomicFree.hypre.arc)
  arcanefem_add_gpu_test(NAME [poisson]2D_bsr_atomicFree_hypre_gpu COMMAND Poisson ARGS inputs/circle.2D.bsr.atomicFree.hypre.arc)

  add_test(NAME [poisson]3D_bsr_hypre COMMAND Poisson inputs/sphere.3D.bsr.hypre.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_bsr_hypre_gpu COMMAND Poisson ARGS inputs/sphere.3D.bsr.hypre.arc)

  add_test(NAME [poisson]3D_neumann_bsr_hypre COMMAND Poisson inputs/sphere.neumann.3D.bsr.hypre.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_neumann_bsr_hypre_gpu COMMAND Poisson ARGS inputs/sphere.neumann.3D.bsr.hypre.arc)

  add_test(NAME [poisson]3D_bsr_atomicFree_hypre COMMAND Poisson inputs/sphere.3D.bsr.atomicFree.hypre.arc)
  arcanefem_add_gpu_test(NAME [poisson]3D_bsr_atomicFree_hypre_gpu COMMAND Poisson ARGS inputs/sphere.3D.bsr.atomicFree.hypre.arc)
endif()

add_test(NAME [poisson]2D_pntDirichlet COMMAND Poisson inputs/perforatedSquare.pointDirichlet.2D.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [poisson]2D_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson inputs/circle.2D.arc)
  add_test(NAME [poisson]3D_4p COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Poisson inputs/sphere.3D.arc)
endif()
