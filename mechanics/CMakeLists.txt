add_executable(Mechanics
  utilFEM.h
  utilFEM.cc
  Types.h
  Mechanics_axl.h
  MechanicsModule.h
  MechanicsModule.cc
  main.cc)

arcane_generate_axl(Mechanics)
arcane_add_arcane_libraries_to_target(Mechanics)

target_include_directories(Mechanics PUBLIC . ../fem1 ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Mechanics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Mechanics PUBLIC FemUtils)
