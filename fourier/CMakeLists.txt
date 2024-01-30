add_executable(Fourier
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Fourier)
target_include_directories(Fourier PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Fourier.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.heterogeneous.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.heterogeneous.10k.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.10k.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.conduction.quad4.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/plancher.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/multi-material.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/plancher.quad4.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Fourier PUBLIC FemUtils)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

find_program(GMSH NAMES gmsh)
if (GMSH)
  message(STATUS "GMSH found: ${GMSH}")
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/plancher.geo -format msh41 -o plancher.10k.msh
    DEPENDS ${MSH_DIR}/plancher.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/multi-material.10k.msh
    COMMAND ${GMSH} -2 ${MSH_DIR}/multi-mat.geo -format msh41 -o multi-material.10k.msh
    DEPENDS ${MSH_DIR}/multi-mat.geo
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )
  add_custom_target(Fourier_gmsh_files DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/multi-material.10k.msh ${CMAKE_CURRENT_BINARY_DIR}/plancher.10k.msh)
  add_dependencies(Fourier Fourier_gmsh_files)
  set(FEMTEST_HAS_GMSH_TEST TRUE)
endif()

enable_testing()

if(FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  add_test(NAME [Fourier]conduction COMMAND Fourier Test.conduction.arc)
  add_test(NAME [Fourier]conduction_heterogeneous COMMAND Fourier Test.conduction.heterogeneous.arc)
  add_test(NAME [Fourier]conduction_quad COMMAND Fourier Test.conduction.quad4.arc)
endif()


# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [Fourier]conduction_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Fourier Test.conduction.arc)
  add_test(NAME [Fourier]conduction_heterogeneous_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Fourier Test.conduction.heterogeneous.arc)
  if(FEMTEST_HAS_GMSH_TEST)
    add_test(NAME [Fourier]conduction_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Fourier Test.conduction.10k.arc)
    add_test(NAME [Fourier]conduction_heterogeneous_10k_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Fourier Test.conduction.heterogeneous.10k.arc)
  endif()
endif()
