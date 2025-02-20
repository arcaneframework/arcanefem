cmake_minimum_required(VERSION 3.21)
project(ArcaneFem LANGUAGES C CXX)

# A partir de 3.27:
# find_package() uses upper-case <PACKAGENAME>_ROOT variables.
if(POLICY CMP0144)
  cmake_policy(SET CMP0144 NEW)
endif()

find_package(Arcane 3.14.14 REQUIRED)
arcane_accelerator_enable()

enable_testing()
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/Utils.cmake)

set(MSH_DIR ${CMAKE_SOURCE_DIR}/meshes/msh)

# Needed for dynamic loading of services
set(BUILD_SHARED_LIBS TRUE)

add_subdirectory(femutils)
add_subdirectory(modules/acoustics)
add_subdirectory(modules/aerodynamics)
add_subdirectory(modules/testlab)
add_subdirectory(modules/laplace)
add_subdirectory(modules/poisson)
add_subdirectory(modules/fourier)
add_subdirectory(modules/electrostatics)
add_subdirectory(modules/bilaplacian)
add_subdirectory(modules/elasticity)
add_subdirectory(modules/elastodynamics)
add_subdirectory(modules/soildynamics)
add_subdirectory(modules/heat)
add_subdirectory(modules/passmo)
add_subdirectory(modules/passmonl)
