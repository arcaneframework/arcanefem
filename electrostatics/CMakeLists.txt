add_executable(Electrostatics
  FemModule.cc
  main.cc
  Fem_axl.h
)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Electrostatics)
target_include_directories(Electrostatics PUBLIC . ${CMAKE_CURRENT_BINARY_DIR})
configure_file(Electrostatics.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Electrostatics.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(Test.Electrostatics.rod-circle.arc ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/box-rods.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
configure_file(${MSH_DIR}/box-rod-circle.msh ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Electrostatics PUBLIC FemUtils)

# Copy the tests files in the binary directory
# The '/' after 'tests' is needed because we want to copy the files
# inside the 'tests' directory but not the directory itself.
file(COPY "tests/" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

enable_testing()

add_test(NAME [electrostatics] COMMAND Electrostatics Test.Electrostatics.arc)
add_test(NAME [electrostatics]rod-circle COMMAND Electrostatics Test.Electrostatics.rod-circle.arc)

# If parallel part is available, add some tests
if(FEMUTILS_HAS_PARALLEL_SOLVER AND MPIEXEC_EXECUTABLE)
  add_test(NAME [electrostatics]parallel_4pe COMMAND ${MPIEXEC_EXECUTABLE} -n 4 ./Electrostatics Test.Electrostatics.arc)
endif()
