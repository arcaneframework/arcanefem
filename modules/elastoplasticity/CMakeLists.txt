set(ACCELERATOR_SOURCES FemModule.cc)

add_executable(Elastoplasticity
        ${ACCELERATOR_SOURCES}
        main.cc
        Fem_axl.h
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Elastoplasticity)

arcane_generate_axl(Fem)
arcane_add_arcane_libraries_to_target(Elastoplasticity)

target_include_directories(Elastoplasticity PUBLIC . ../fem ${CMAKE_CURRENT_BINARY_DIR})
target_include_directories(Elastoplasticity PUBLIC /home/catC/aa286863/exaquake/psd/src/PSDinstall/include)

configure_file(Elastoplasticity.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elastoplasticity PUBLIC FemUtils /home/catC/aa286863/exaquake/psd/src/PSDinstall/lib/libMFrontGenericInterface.so)

#  Copy the data files in the binary directory
file(COPY "data" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the check files in the binary directory
file(COPY "check" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the inputs folder containing the arc files
file(COPY "inputs" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# Copy the mesh files into meshes directory
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/meshes)
set(MESH_FILES
        quater_cylinder.msh
        quater_cylinder-quad4.msh
        quater_cylinder-quad8.msh
        quater_cylinder-quad9.msh
        quater_cylinder-multi_domain.msh
        strip_footing.msh
)

foreach (MESH_FILE IN LISTS MESH_FILES)
  file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach ()

enable_testing()
set(MODULE_NAME "elastoplasticity")

if (FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
	set(SOLVER_PETSC_GMRES
		-A,//fem/petsc-flags=-ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15\ -ksp_atol\ 1e-15)

  set(SOLVER_ALEPH_CG
		-A,//fem/linear-system/@name=AlephLinearSystem
		-A,//fem/linear-system/solver=cg
		-A,//fem/linear-system/epsilon=1.0e-12)

  set(SOLVER_ALEPH_CG_GAMG
    -A,//fem/linear-system/@name=AlephLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/preconditioner=amg
    -A,//fem/linear-system/amg-near-null-space=true
    -A,//fem/petsc-flags=-ksp_type\ cg\ -mat_block_size\ 2\ -ksp_rtol\ 1e-8\ -pc_type\ gamg\ -pc_gamg_threshold\ 0.001
    -A,//fem/linear-system/epsilon=1.0e-8)

  arcanefem_add_full_test(NAME petsc_gamg_near_null_space
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//fem/matrix-format=BSR
      -A,//fem/linear-system/@name=PetscLinearSystem
      -A,//fem/linear-system/solver=cg
      -A,//fem/linear-system/pc-type=gamg
      -A,//fem/linear-system/amg-near-null-space=true
      -A,//fem/linear-system/rtol=1.0e-9
      -A,//fem/petsc-flags=-mat_block_size\ 2
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      ${SOLVER_ALEPH_CG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_parse_and_exit
    COMMAND ./Elastoplasticity
    ARGS
		  -A,//fem/solve-nonlinear-system=false
		  inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_two_domains_von_mises_dok_global
    COMMAND ./Elastoplasticity
    ARGS
      ${SOLVER_ALEPH_CG}
      inputs/2D.quater.cylinder.two-domains.vonMises.arc)

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad4
    COMMAND ./Elastoplasticity
    ARGS
      -A,//meshes/mesh/filename=meshes/quater_cylinder-quad4.msh
      -A,//fem/matrix-format=DOK
      -A,//fem/hex-quad-mesh=true
      ${SOLVER_ALEPH_CG_GAMG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad8
    COMMAND ./Elastoplasticity
    ARGS
      -A,//meshes/mesh/filename=meshes/quater_cylinder-quad8.msh
      -A,//fem/matrix-format=DOK
      -A,//fem/hex-quad-mesh=true
      ${SOLVER_ALEPH_CG_GAMG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad9
    COMMAND ./Elastoplasticity
    ARGS
      -A,//meshes/mesh/filename=meshes/quater_cylinder-quad9.msh
      -A,//fem/matrix-format=DOK
      -A,//fem/hex-quad-mesh=true
      ${SOLVER_ALEPH_CG_GAMG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_ALEPH_CG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr
    COMMAND ./Elastoplasticity
    ARGS
		  -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
		  -A,//fem/matrix-format=AF-BSR
		  ${SOLVER_ALEPH_CG}
		  inputs/2D.quater.cylinder.vonMises.arc)

  set(SOLVER_PETSC_DIRECT_CG
    -A,//fem/linear-system/@name=PetscLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/max-iter=5000
    -A,//fem/linear-system/atol=1e-12
    -A,//fem/linear-system/rtol=0
    -A,//fem/petsc-flags=-ksp_view\ -mat_block_size\ 2\ -pc_type\ gamg
    )

  set(SOLVER_PETSC_DIRECT_GMRES
    -A,//fem/linear-system/@name=PetscLinearSystem
    -A,//fem/linear-system/solver=gmres
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0
    -A,//fem/petsc-flags=-ksp_view\ -mat_block_size\ 2\ -pc_type\ gamg
    )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_bsr_petsc-direct
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_DIRECT_CG}
      inputs/2D.strip.footing.DruckerPrager.arc)

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_af-bsr_petsc-direct
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_DIRECT_CG}
      inputs/2D.strip.footing.DruckerPrager.arc)

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr_petsc-direct
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_PETSC_DIRECT_CG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr_petsc-direct
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_PETSC_DIRECT_CG}
      inputs/2D.quater.cylinder.vonMises.arc)
endif ()

if (FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=cg
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)

  set(SOLVER_HYPRE_BICGSTAB
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=bicgstab)

  set(SOLVER_HYPRE_GMRES
    -A,//fem/linear-system/@name=HypreLinearSystem
    -A,//fem/linear-system/solver=gmres
    -A,//fem/linear-system/atol=1e-9
    -A,//fem/linear-system/rtol=0)

  arcanefem_add_seq_test(NAME hypre_amg_near_null_space
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//fem/matrix-format=BSR
      -A,//fem/linear-system/@name=HypreLinearSystem
      -A,//fem/linear-system/solver=cg
      -A,//fem/linear-system/preconditioner=amg
      -A,//fem/linear-system/amg-near-null-space=true
      -A,//fem/linear-system/rtol=1.0e-9
      inputs/2D.quater.cylinder.vonMises.arc)

  if (FEMUTILS_HAS_PARALLEL_SOLVER)
    arcanefem_add_mpi_test(NAME hypre_amg_near_null_space
      COMMAND ./Elastoplasticity
      ARGS
        -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
        -A,//fem/matrix-format=BSR
        -A,//fem/linear-system/@name=HypreLinearSystem
        -A,//fem/linear-system/solver=cg
        -A,//fem/linear-system/preconditioner=amg
        -A,//fem/linear-system/amg-near-null-space=true
        -A,//fem/linear-system/rtol=1.0e-9
        inputs/2D.quater.cylinder.vonMises.arc)
  endif ()

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr_hypre
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr_hypre
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG}
      inputs/2D.quater.cylinder.vonMises.arc)

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_bsr_hypre
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
      -A,//fem/matrix-format=BSR
      ${SOLVER_HYPRE_CG}
      inputs/2D.strip.footing.DruckerPrager.arc)

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_af-bsr_hypre
    COMMAND ./Elastoplasticity
    ARGS
      -A,//fem/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
      -A,//fem/matrix-format=AF-BSR
      ${SOLVER_HYPRE_CG}
      inputs/2D.strip.footing.DruckerPrager.arc)
endif ()
