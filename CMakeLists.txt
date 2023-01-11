cmake_minimum_required(VERSION 3.21)
project(FemTest1Main LANGUAGES NONE)
enable_testing()

set(MSH_DIR ${CMAKE_SOURCE_DIR}/meshes/msh)

add_subdirectory(femutils)
add_subdirectory(poisson)
add_subdirectory(bilaplacian)
add_subdirectory(elasticity)
add_subdirectory(elastodynamics)
add_subdirectory(heat)
add_subdirectory(passmo)
