set(ACCELERATOR_SOURCES
  Elastoplasticity2Module.cc
  Dirichlet.cc
  DruckerPragerLaw.cc
  VonMisesLaw.cc
  InternalBodyForce.cc
)

add_executable(Elastoplasticity2
  ${ACCELERATOR_SOURCES}
  Elastoplasticity2Module.h
  Elastoplasticity2_axl.h
  ElementMatrix.h
  ElementMatrixHexQuad.h
  ExternalBodyForce.cc
  Traction.cc
  main.cc
)

arcane_accelerator_add_source_files(${ACCELERATOR_SOURCES})
arcane_accelerator_add_to_target(Elastoplasticity2)

arcane_generate_axl(Elastoplasticity2)
arcane_add_arcane_libraries_to_target(Elastoplasticity2)

target_include_directories(Elastoplasticity2 PUBLIC . ${CMAKE_BINARY_DIR} ${ArcaneFem_SOURCE_DIR})

configure_file(Elastoplasticity2.config ${CMAKE_CURRENT_BINARY_DIR} COPYONLY)

target_link_libraries(Elastoplasticity2 PUBLIC FemUtils)

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
  strip_footing.msh
)

foreach (MESH_FILE IN LISTS MESH_FILES)
  file(COPY ${MSH_DIR}/${MESH_FILE} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/meshes)
endforeach ()

enable_testing()
set(MODULE_NAME "elastoplasticity2")

if (FEMUTILS_HAS_SOLVER_BACKEND_PETSC)
  set(SOLVER_PETSC_GMRES
    -A,//elastoplasticity2/petsc-flags=-ksp_monitor\ -ksp_type\ gmres\ -pc_type\ bjacobi\ -ksp_rtol\ 1e-15\ -ksp_atol\ 1e-15)

  set(SOLVER_ALEPH_CG
    -A,//elastoplasticity2/linear-system/@name=AlephLinearSystem
    -A,//elastoplasticity2/linear-system/solver=cg
    -A,//elastoplasticity2/linear-system/epsilon=1.0e-12)

  set(SOLVER_ALEPH_CG_GAMG
    -A,//elastoplasticity2/linear-system/@name=AlephLinearSystem
    -A,//elastoplasticity2/linear-system/solver=cg
    -A,//elastoplasticity2/linear-system/preconditioner=amg
    -A,//elastoplasticity2/linear-system/amg-near-null-space=true
    -A,//elastoplasticity2/petsc-flags=-ksp_type\ cg\ -mat_block_size\ 2\ -ksp_rtol\ 1e-8\ -pc_type\ gamg\ -pc_gamg_threshold\ 0.001
    -A,//elastoplasticity2/linear-system/epsilon=1.0e-8)

  arcanefem_add_full_test(NAME petsc_gamg_near_null_space
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    -A,//elastoplasticity2/linear-system/@name=PetscLinearSystem
    -A,//elastoplasticity2/linear-system/solver=cg
    -A,//elastoplasticity2/linear-system/pc-type=gamg
    -A,//elastoplasticity2/linear-system/amg-near-null-space=true
    -A,//elastoplasticity2/linear-system/rtol=1.0e-9
    -A,//elastoplasticity2/petsc-flags=-mat_block_size\ 2
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    ${SOLVER_ALEPH_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_local
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_ALEPH_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_parse_and_exit
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solve-nonlinear-system=false
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad4
    COMMAND ./Elastoplasticity2 ARGS
    -A,//meshes/mesh/filename=meshes/quater_cylinder-quad4.msh
    -A,//elastoplasticity2/matrix-format=DOK
    -A,//elastoplasticity2/hex-quad-mesh=true
    ${SOLVER_ALEPH_CG_GAMG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad4_local
    COMMAND ./Elastoplasticity2 ARGS
    -A,//meshes/mesh/filename=meshes/quater_cylinder-quad4.msh
    -A,//elastoplasticity2/matrix-format=DOK
    -A,//elastoplasticity2/hex-quad-mesh=true
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_ALEPH_CG_GAMG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad8
    COMMAND ./Elastoplasticity2 ARGS
    -A,//meshes/mesh/filename=meshes/quater_cylinder-quad8.msh
    -A,//elastoplasticity2/matrix-format=DOK
    -A,//elastoplasticity2/hex-quad-mesh=true
    ${SOLVER_ALEPH_CG_GAMG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad8_local
    COMMAND ./Elastoplasticity2 ARGS
    -A,//meshes/mesh/filename=meshes/quater_cylinder-quad8.msh
    -A,//elastoplasticity2/matrix-format=DOK
    -A,//elastoplasticity2/hex-quad-mesh=true
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_ALEPH_CG_GAMG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad9
    COMMAND ./Elastoplasticity2 ARGS
    -A,//meshes/mesh/filename=meshes/quater_cylinder-quad9.msh
    -A,//elastoplasticity2/matrix-format=DOK
    -A,//elastoplasticity2/hex-quad-mesh=true
    ${SOLVER_ALEPH_CG_GAMG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_seq_test(NAME 2D_Quater_cylinder_von_mises_quad9_local
    COMMAND ./Elastoplasticity2 ARGS
    -A,//meshes/mesh/filename=meshes/quater_cylinder-quad9.msh
    -A,//elastoplasticity2/matrix-format=DOK
    -A,//elastoplasticity2/hex-quad-mesh=true
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_ALEPH_CG_GAMG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    ${SOLVER_ALEPH_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    ${SOLVER_ALEPH_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr_local
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_ALEPH_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr_local
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_ALEPH_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  set(SOLVER_PETSC_DIRECT_CG
    -A,//elastoplasticity2/linear-system/@name=PetscLinearSystem
    -A,//elastoplasticity2/linear-system/solver=cg
    -A,//elastoplasticity2/linear-system/max-iter=5000
    -A,//elastoplasticity2/linear-system/atol=1e-12
    -A,//elastoplasticity2/linear-system/rtol=0
    -A,//elastoplasticity2/petsc-flags=-ksp_view\ -mat_block_size\ 2\ -pc_type\ gamg
  )

  set(SOLVER_PETSC_DIRECT_GMRES
    -A,//elastoplasticity2/linear-system/@name=PetscLinearSystem
    -A,//elastoplasticity2/linear-system/solver=gmres
    -A,//elastoplasticity2/linear-system/atol=1e-9
    -A,//elastoplasticity2/linear-system/rtol=0
    -A,//elastoplasticity2/petsc-flags=-ksp_view\ -mat_block_size\ 2\ -pc_type\ gamg
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_bsr_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=BSR
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_af-bsr_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_bsr_local_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_af-bsr_local_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr_local_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr_local_petsc-direct
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_PETSC_DIRECT_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )
endif()

if (FEMUTILS_HAS_SOLVER_BACKEND_HYPRE)
  set(SOLVER_HYPRE_CG
    -A,//elastoplasticity2/linear-system/@name=HypreLinearSystem
    -A,//elastoplasticity2/linear-system/solver=cg
    -A,//elastoplasticity2/linear-system/atol=1e-9
    -A,//elastoplasticity2/linear-system/rtol=0)

  set(SOLVER_HYPRE_BICGSTAB
    -A,//elastoplasticity2/linear-system/@name=HypreLinearSystem
    -A,//elastoplasticity2/linear-system/solver=bicgstab)

  set(SOLVER_HYPRE_GMRES
    -A,//elastoplasticity2/linear-system/@name=HypreLinearSystem
    -A,//elastoplasticity2/linear-system/solver=gmres
    -A,//elastoplasticity2/linear-system/atol=1e-9
    -A,//elastoplasticity2/linear-system/rtol=0)

  arcanefem_add_seq_test(NAME hypre_amg_near_null_space
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    -A,//elastoplasticity2/linear-system/@name=HypreLinearSystem
    -A,//elastoplasticity2/linear-system/solver=cg
    -A,//elastoplasticity2/linear-system/preconditioner=amg
    -A,//elastoplasticity2/linear-system/amg-near-null-space=true
    -A,//elastoplasticity2/linear-system/rtol=1.0e-9
    inputs/2D.quater.cylinder.vonMises.arc
  )

  if (FEMUTILS_HAS_PARALLEL_SOLVER)
    arcanefem_add_mpi_test(NAME hypre_amg_near_null_space
      COMMAND ./Elastoplasticity2 ARGS
      -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
      -A,//elastoplasticity2/matrix-format=BSR
      -A,//elastoplasticity2/linear-system/@name=HypreLinearSystem
      -A,//elastoplasticity2/linear-system/solver=cg
      -A,//elastoplasticity2/linear-system/preconditioner=amg
      -A,//elastoplasticity2/linear-system/amg-near-null-space=true
      -A,//elastoplasticity2/linear-system/rtol=1.0e-9
      inputs/2D.quater.cylinder.vonMises.arc
    )
  endif ()

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr_hypre
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr_hypre
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_bsr_local_hypre
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_HYPRE_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Quater_cylinder_von_mises_af-bsr_local_hypre
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.quater.cylinder.t20.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_HYPRE_CG}
    inputs/2D.quater.cylinder.vonMises.arc
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_bsr_hypre
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=BSR
    ${SOLVER_HYPRE_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_af-bsr_hypre
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    ${SOLVER_HYPRE_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_bsr_local_hypre
    COMMAND ./Elastoplasticity2  ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_HYPRE_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )

  arcanefem_add_full_test(NAME 2D_Strip_footing_drucker_prager_af-bsr_local_hypre
    COMMAND ./Elastoplasticity2 ARGS
    -A,//elastoplasticity2/solution-comparison-file=check/2D.strip.footing.drucker_prager.t12.txt
    -A,//elastoplasticity2/matrix-format=AF-BSR
    -A,//elastoplasticity2/gp-material-tensor-strategy=local
    ${SOLVER_HYPRE_CG}
    inputs/2D.strip.footing.DruckerPrager.arc
  )
endif()

# ----------------------------------------------------------------------------
# Local Variables:
# tab-width: 2
# indent-tabs-mode: nil
# coding: utf-8-with-signature
# End:
