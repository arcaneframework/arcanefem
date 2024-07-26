add_executable(Acoustics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Acoustics)
target_include_directories(Acoustics PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Acoustics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.sub.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/sub.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Acoustics PUBLIC FemUtils)

enable_testing()

add_test(NAME [Acoustics]2D_submarine COMMAND Acoustics Test.sub.arc)
