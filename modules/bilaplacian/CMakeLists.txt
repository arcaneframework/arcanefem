cmake_minimum_required(VERSION 3.21)
project(Bilaplacian LANGUAGES C CXX)

find_package(Arcane REQUIRED)

add_executable(Bilaplacian
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Bilaplacian)
target_include_directories(Bilaplacian PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Bilaplacian.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Bilaplacian PUBLIC FemUtils)

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  bilap.msh
)
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

add_test(NAME [bilaplacian]direct_solver COMMAND Bilaplacian inputs/direct.arc)
add_test(NAME [bilaplacian]parse_and_exit COMMAND Bilaplacian
  -A,//fem/assemble-linear-system=false
  -A,//fem/solve-linear-system=false
  -A,//fem/cross-validation=false
  inputs/direct.arc)

#if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
#  add_test(NAME [bilaplacian]petsc COMMAND Bilaplacian
#    -A,//fem/linear-system/solver-backend=petsc
#    -A,//fem/linear-system/solver-backend=petsc -A,//fem/petsc-flags=-ksp_monitor\ -ksp_type\ gmres\ -pc_type\ lu\ -ksp_rtol\ 1e-15\ -mat_type\ mpibaij\ -bs\ 2
#    inputs/internal_hypre_rowColElim.arc)
#endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [bilaplacian]hypre_rowColElim COMMAND Bilaplacian
    inputs/internal_hypre_rowColElim.arc)
  add_test(NAME [bilaplacian]hypre_rowElim COMMAND Bilaplacian
    -A,//fem/enforce-Dirichlet-method=RowElimination
    inputs/internal_hypre_rowColElim.arc)
endif()
