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
add_subdirectory(acoustics)
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
