set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Fourier
  FemModule.h
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Fourier)
target_include_directories(Fourier PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Fourier.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Fourier PUBLIC FemUtils)


if (NOT APPLE)
  # Command to compile C# file for functions
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ExternalFunctions.dll"
    COMMAND ${ARCANE_PREFIX_DIR}/bin/arcane_dotnet_compile "${CMAKE_CURRENT_LIST_DIR}/ExternalFunctions.cs"
    DEPENDS "${CMAKE_CURRENT_LIST_DIR}/ExternalFunctions.cs"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  add_custom_target(external_functions DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/ExternalFunctions.dll")
  add_dependencies(Fourier external_functions)
else()
  message(STATUS "Skipping ExternalFunctions.dll build on macOS")
endif()

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  plancher.msh
  plancher.quad4.msh
  multi-material.msh
  multi-material.quad.msh
  square_-2pi_to_2pi.msh
  sphere_cut.hexa.msh
  bar_dynamic_3D.msh
  pipe.msh
  circle_cut.quad.msh
  boolean.msh
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
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/meshes/multi-material.10k.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/multi-mat.geo -format msh41 -o meshes/multi-material.10k.msh -bin
    DEPENDS ${MSH_DIR}/multi-mat.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(Fourier_gmsh_files DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/meshes/multi-material.10k.msh ${CMAKE_CURRENT_BINARY_DIR}/meshes/plancher.10k.msh)
  add_dependencies(Fourier Fourier_gmsh_files)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [fourier]conduction COMMAND Fourier
    -A,//fem/solution-comparison-file=check/test1_results.txt
    inputs/conduction.arc)
  add_test(NAME [fourier]conduction_heterogeneous COMMAND Fourier
    -A,//fem/solution-comparison-file=check/test2_results.txt
    inputs/conduction.heterogeneous.arc)
  add_test(NAME [fourier]conduction_quad COMMAND Fourier
    -A,//fem/solution-comparison-file=check/conduction_quad.txt
    inputs/conduction.quad.arc)
  add_test(NAME [fourier]conduction_heterogeneous_quad COMMAND Fourier
    -A,//fem/solution-comparison-file=check/conduction_multi-mat_quad.txt
    inputs/conduction.heterogeneous.quad.arc)
  add_test(NAME [fourier]conduction_3D COMMAND Fourier
    -A,//fem/solution-comparison-file=check/test_conduction_3D.txt
    inputs/conduction.3D.arc)
  add_test(NAME [fourier]conduction_boolean_3D COMMAND Fourier
    inputs/conduction.boolean.3D.arc)
  add_test(NAME [fourier]conduction_sphere-cut_hexa COMMAND Fourier
    -A,//fem/solution-comparison-file=check/conduction_hexa.txt
    inputs/conduction.sphere-cut.hexa.arc)

  if (NOT APPLE)
    add_test(NAME [fourier]manufacture_solution COMMAND Fourier
      -A,UsingDotNet=1
      inputs/manufacture.solution.arc)
    add_test(NAME [fourier]manufacture_solution_3D COMMAND Fourier
      -A,UsingDotNet=1
      -A,//meshes/mesh/filename=meshes/bar_dynamic_3D.msh
      inputs/manufacture.solution.arc)
  endif()

  arcanefem_add_gpu_test(NAME [fourier]conduction_petsc_bsr COMMAND ./Fourier ARGS
    -A,//fem/solution-comparison-file=check/test1_results.txt
    -A,//fem/matrix-format=BSR
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [fourier]conduction_petsc_af-bsr COMMAND ./Fourier ARGS
    -A,//fem/solution-comparison-file=check/test1_results.txt
    -A,//fem/matrix-format=AF-BSR
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [fourier]3d_conduction_petsc_bsr COMMAND ./Fourier ARGS
    -A,//fem/matrix-format=BSR
    inputs/conduction.3D.arc)

  arcanefem_add_gpu_test(NAME [fourier]3d_conduction_petsc_af-bsr COMMAND ./Fourier ARGS
    -A,//fem/matrix-format=AF-BSR
    inputs/conduction.3D.arc)

  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [fourier]conduction_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
      -A,//fem/solution-comparison-file=check/test1_results.txt
      inputs/conduction.arc)
    add_test(NAME [fourier]conduction_heterogeneous_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
      -A,//fem/solution-comparison-file=check/test2_results.txt
      inputs/conduction.heterogeneous.arc)
    add_test(NAME [fourier]conduction_quad_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
      -A,//fem/solution-comparison-file=check/conduction_quad.txt
      inputs/conduction.quad.arc)
    add_test(NAME [fourier]conduction_heterogeneous_quad_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
      -A,//fem/solution-comparison-file=check/conduction_multi-mat_quad.txt
      inputs/conduction.heterogeneous.quad.arc)
    add_test(NAME [fourier]conduction_heterogeneous_3D_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
      inputs/conduction.pipe.3D.arc)
    add_test(NAME [fourier]conduction_sphere-cut_3D_hexa_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
      -A,//fem/solution-comparison-file=check/conduction_hexa.txt
      inputs/conduction.sphere-cut.hexa.arc)

    if (NOT APPLE)
      add_test(NAME [fourier]manufacture_solution_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
        -A,UsingDotNet=1
        inputs/manufacture.solution.arc)
      add_test(NAME [fourier]manufacture_solution_3D_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
        -A,UsingDotNet=1
        -A,//meshes/mesh/filename=meshes/bar_dynamic_3D.msh
        inputs/manufacture.solution.arc)
    endif()

    if(FEMTEST_HAS_GMSH_TEST)
      add_test(NAME [fourier]conduction_10k_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
        -A,//meshes/mesh/filename=meshes/plancher.10k.msh
        -A,//fem/cross-validation=false
        inputs/conduction.arc)
      add_test(NAME [fourier]conduction_heterogeneous_10k_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Fourier
        -A,//meshes/mesh/filename=meshes//multi-material.10k.msh
        -A,//fem/cross-validation=false
        inputs/conduction.heterogeneous.arc)
    endif()
  endif()
endif()

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_gpu_test(NAME [fourier]conduction_hypre_bsr COMMAND ./Fourier ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [fourier]conduction_hypre_af-bsr COMMAND ./Fourier ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/conduction.arc)

  arcanefem_add_gpu_test(NAME [fourier]3d_conduction_hypre_bsr COMMAND ./Fourier ARGS
    -A,//fem/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/conduction.3D.arc)

  arcanefem_add_gpu_test(NAME [fourier]3d_conduction_hypre_af-bsr COMMAND ./Fourier ARGS
    -A,//fem/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/conduction.3D.arc)

endif()