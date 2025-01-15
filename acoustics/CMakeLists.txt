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
configure_file(Test.sub.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.sub.hypre.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/sub.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(Acoustics PUBLIC FemUtils)

enable_testing()

add_test(NAME [Acoustics]2D_submarine COMMAND Acoustics Test.sub.arc)

if(FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  add_test(NAME [Acoustics]2D_submarine_hypre COMMAND Acoustics Test.sub.hypre.arc)
  if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
    add_test(NAME [Acoustics]2D_submarine_hypre_2pe COMMAND ${MPIEXEC_EXECUTABLE} -n 2 ./Acoustics Test.sub.hypre.arc)
  endif()
endif()
