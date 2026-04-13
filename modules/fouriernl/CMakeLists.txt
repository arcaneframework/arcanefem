set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(FourierNL
  FemModule.h
  ${ACCELERATOR_SOURCES}
  main.cc
  Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(FourierNL)
target_include_directories(FourierNL PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(FourierNL.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(FourierNL PUBLIC FemUtils)


find_program(ARCANE_DOTNET_COMPILE
  NAMES arcane_dotnet_compile
  PATHS ${ARCANE_PREFIX_DIR}/bin
)

if (NOT APPLE AND ARCANE_DOTNET_COMPILE
    AND EXISTS "${ARCANE_PREFIX_DIR}/lib/Arcane.Core.dll"
    AND EXISTS "${ARCANE_PREFIX_DIR}/lib/Arcane.Utils.dll")
  message(STATUS "Arcane .NET support detected — building ExternalFunctions.dll")

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/ExternalFunctions.dll"
    COMMAND ${ARCANE_DOTNET_COMPILE} "${CMAKE_CURRENT_LIST_DIR}/ExternalFunctions.cs"
    DEPENDS "${CMAKE_CURRENT_LIST_DIR}/ExternalFunctions.cs"
    WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  )

  add_custom_target(external_functions DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/ExternalFunctions.dll")
  add_dependencies(FourierNL external_functions)

elseif(APPLE)
  message(STATUS "Skipping ExternalFunctions.dll build on macOS")

else()
  message(STATUS "Skipping ExternalFunctions.dll — Arcane built without .NET support")
endif()

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
  square_-2pi_to_2pi.msh
  )
foreach(MESH_FILE IN LISTS MESH_FILES)
    file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  if (NOT APPLE)
    add_test(NAME [fouriernl]manufacture_solution COMMAND FourierNL
      -A,UsingDotNet=1
      inputs/manufacture.solution.arc)
    add_test(NAME [fouriernl]manufacture_solution_3D COMMAND FourierNL
      -A,UsingDotNet=1
      -A,//meshes/mesh/filename=meshes/bar_dynamic_3D.msh
      inputs/manufacture.solution.arc)
  endif()

  set(SOLVER_PETSC_CG
    -A,//fem/linear-system/@name=PETScLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)


  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    if (NOT APPLE)
      add_test(NAME [fouriernl]manufacture_solution_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./FourierNL
        -A,UsingDotNet=1
        inputs/manufacture.solution.arc)
      add_test(NAME [fouriernl]manufacture_solution_3D_2p COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./FourierNL
        -A,UsingDotNet=1
        -A,//meshes/mesh/filename=meshes/bar_dynamic_3D.msh
        inputs/manufacture.solution.arc)
    endif()
  endif()
endif()

#if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
#  set(SOLVER_HYPRE_CG
#    -A,//fem/linear-system/@name=HypreLinearSystem
#    -A,//fem/linear-system/solver=cg
#    -A,//fem/linear-system/atol=1e-9
#    -A,//fem/linear-system/rtol=0)
#
#  arcanefem_add_gpu_test(NAME [fourier]conduction_hypre_bsr COMMAND ./Fourier ARGS
#    -A,//fem/matrix-format=BSR
#    ${SOLVER_HYPRE_CG}
#    inputs/conduction.arc)
#
#  arcanefem_add_gpu_test(NAME [fourier]conduction_hypre_af-bsr COMMAND ./Fourier ARGS
#    -A,//fem/matrix-format=AF-BSR
#    ${SOLVER_HYPRE_CG}
#    inputs/conduction.arc)
#
#  arcanefem_add_gpu_test(NAME [fourier]3d_conduction_hypre_bsr COMMAND ./Fourier ARGS
#    -A,//fem/matrix-format=BSR
#    ${SOLVER_HYPRE_CG}
#    inputs/conduction.3D.arc)
#
#  arcanefem_add_gpu_test(NAME [fourier]3d_conduction_hypre_af-bsr COMMAND ./Fourier ARGS
#    -A,//fem/matrix-format=AF-BSR
#    ${SOLVER_HYPRE_CG}
#    inputs/conduction.3D.arc)
#
#endif()