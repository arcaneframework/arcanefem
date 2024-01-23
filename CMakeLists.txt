cmake_minimum_required(VERSION 3.21)
project(ArcaneFem LANGUAGES C CXX)

find_package(Arcane REQUIRED)

enable_testing()

set(MSH_DIR ${CMAKE_SOURCE_DIR}/meshes/msh)

# Needed for dynamic loading of services
set(BUILD_SHARED_LIBS TRUE)

add_subdirectory(femutils)
add_subdirectory(aerodynamics)
add_subdirectory(poisson)
add_subdirectory(laplace)
add_subdirectory(fourier)
add_subdirectory(electrostatics)
add_subdirectory(bilaplacian)
add_subdirectory(elasticity)
add_subdirectory(elastodynamics)
add_subdirectory(soildynamics)
add_subdirectory(heat)
add_subdirectory(passmo)
